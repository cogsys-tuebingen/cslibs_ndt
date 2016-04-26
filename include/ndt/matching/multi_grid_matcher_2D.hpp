#ifndef MATCHER2D_HPP
#define MATCHER2D_HPP

#include <ndt/grid/multi_grid.hpp>
#include <ndt/matching/matcher.hpp>
#include <ndt/data/pointcloud.hpp>

#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace ndt {
namespace matching {
class MultiGridMatcher2D : public Matcher<2> {
public:
    typedef grid::MultiGrid<2>                         GridType;
    typedef typename GridType::DistributionType        DistributionType;
    typedef typename GridType::DistributionSetType     DistributionsType;
    typedef typename GridType::CovarianceMatrixType    CovarianceMatrixType;
    typedef typename GridType::SizeType                SizeType;
    typedef typename PointCloudType::PointType         PointType;
    typedef Eigen::Matrix<double,3,3>                  HessianType;
    typedef Eigen::Translation2d                       TranslationType;
    typedef Eigen::Rotation2Dd                         RotationType;
    typedef Eigen::Vector3d                            GradientType;


    MultiGridMatcher2D(const Parameters &params = Parameters()) :
        BaseClass(params)
    {
    }

    virtual ~MultiGridMatcher2D()
    {
    }

    inline double match(const PointCloudType &_src,
                        const PointCloudType &_dst,
                        TransformType        &_transformation,
                        const TransformType  &_prior_transformation = TransformType::Identity()) override
    {
        /// build the ndt grid for the src cloud
        PointType range = _src.range();
        typename GridType::SizeType size;
        for(std::size_t i = 0 ; i < 2 ; ++i) {
            if(range(i) <= 0.0)
                throw std::range_error("Point cloud boundaries are not set properly!");
            size[i] = floor(range(i) / params.resolution[i] + 0.5);
        }

        grid.reset(new GridType(size, params.resolution, _src.min));
        grid->add(_src);

        /// todo:: initialize parameter estimate double phi = 0.0;
        double          tx = _prior_transformation.translation()(0);
        double          ty = _prior_transformation.translation()(1);
        double          phi = acos(_prior_transformation.rotation()(0,0));

        TranslationType trans;
        RotationType    rotation(0.0);

        /// variables needed for sampling a point
        double          sin_phi;
        double          cos_phi;

        /// gradient and stuff
        /// need 4 fields for that
        /// only solve the maximal score
        double                      max_score;
        std::size_t                 max_idx;
        std::array<double, 4>       score;
        std::array<GradientType, 4> gradient;
        std::array<HessianType, 4>  hessian;
        GradientType                delta_p;

        double       tx_old;
        double       ty_old;
        double       phi_old;
        std::size_t  iteration = 0;

        bool converged = false;
        while(!converged) {
            rotation = RotationType(phi);
            trans    = TranslationType(tx, ty);
            _transformation = trans * rotation;

            gradient.fill(GradientType::Zero());
            hessian.fill(HessianType::Zero());
            score.fill(0.0);
            tx_old   = tx;
            ty_old   = ty;
            phi_old  = phi;
            sincos(phi, &sin_phi, &cos_phi);

            /// calculate the hessian and the gradient for each grid
            for(std::size_t i = 0 ; i < _dst.size ; ++i) {
                if(_dst.mask[i] == PointCloudType::VALID) {
                    PointType p = _transformation * _dst.points[i];

                    DistributionsType distributions;
                    grid->get(p, distributions);
                    for(std::size_t j = 0 ; j < distributions.size(); ++j) {
                        if(distributions[j] == nullptr)
                            continue;
                        DistributionType &distribution = *distributions[j];
                        if(distribution.getN() < 3)
                            continue;

                        PointType            mean;
                        PointType            q;
                        CovarianceMatrixType inverse_covariance;
                        PointType            jac;
                        PointType            hes;
                        GradientType        &gradient_entry = gradient[j];
                        HessianType         &hessian_entry  = hessian[j];

                        double s = distribution.sampleNonNoramlized(p, q);
                        distribution.getMean(mean);
                        distribution.getInverseCovariance(inverse_covariance);

                        /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                        /// if s is 0.0 we do no need to count the sample in
                        if(s > 1e-3) {
                            score[j] += s;
                            const PointType q_inverse_covariance = q.transpose() * inverse_covariance;

                            jac(0) = -q(0) * sin_phi - q(1) * cos_phi;
                            jac(1) =  q(0) * cos_phi - q(1) * sin_phi;
                            hes(0) = -q(0) * cos_phi + q(1) * sin_phi;
                            hes(1) = -q(0) * sin_phi - q(1) * cos_phi;

                            /// gradient computation
                            double g_dot = q_inverse_covariance.dot(jac);
                            gradient_entry(0) -= s * q_inverse_covariance(0);
                            gradient_entry(1) -= s * q_inverse_covariance(1);
                            gradient_entry(2) -= s * g_dot;
                            /// hessian computation
                            hessian_entry(0,0)+=  s
                                    * (-(q_inverse_covariance(0) * q_inverse_covariance(0))   /// (1)
                                       +  inverse_covariance(0,0));                           /// (3)
                            hessian_entry(1,0)+=  s
                                    * (-(q_inverse_covariance(1) * q_inverse_covariance(0))   /// (1)
                                       +  inverse_covariance(1,0));                           /// (3)
                            hessian_entry(2,0)+=  s
                                    * (-(g_dot * q_inverse_covariance(0))                     /// (1)
                                       +(inverse_covariance.row(0).dot(jac)));                /// (3)
                            hessian_entry(0,1)+=  s
                                    * (-(q_inverse_covariance(0) * q_inverse_covariance(1))   /// (1)
                                       +  inverse_covariance(0,1));                           /// (3)
                            hessian_entry(1,1)+=  s
                                    * (-(q_inverse_covariance(1) * q_inverse_covariance(1))   /// (1)
                                       +  inverse_covariance(1,1));                           /// (3)
                            hessian_entry(2,1)+=  s
                                    * (-(g_dot * q_inverse_covariance(1))                     /// (1)
                                       +(inverse_covariance.row(1).dot(jac)));                /// (3)
                            hessian_entry(0,2)+=  s
                                    * (-(q_inverse_covariance(0) * g_dot)                     /// (1)
                                       +(jac.transpose().dot(inverse_covariance.col(0))));    /// (3)
                            hessian_entry(1,2)+=  s
                                    * (-(q_inverse_covariance(1) * g_dot )                    /// (1)
                                       +(jac.transpose() * inverse_covariance.col(1)));       /// (3)
                            hessian_entry(2,2)+=  s
                                    * (-g_dot * g_dot                                         /// (1)
                                       +q_inverse_covariance.dot(hes)                         /// (2)
                                       +(jac.transpose() * inverse_covariance * jac));        /// (3)

                            /// (1) directly computed from Jacobian combined with q^t * InvCov
                            /// (2) only a result for H(2,2)
                            /// (3) [1,0].[[a,b],[c,d]].[[1],[0]] = a with i = j = 0, Jac.col(0)
                            ///     [0,1].[[a,b],[c,d]].[[1],[0]] = c with i = 1, j = 0, Jac.col(1), Jac.col(0)
                            ///     [1,0].[[a,b],[c,d]].[[0],[1]] = b with i = 0, j = 1, Jac.col(0), Jac.col(1)
                            ///     [0,1].[[a,b],[c,d]].[[0],[1]] = d with i = 1, j = 1, Jac.col(1)
                            ///     [1,0].[[a,b],[c,d]].J_T.col(2) = [a,b].J_T.col(2)
                            ///     [0,1].[[a,b],[c,d]].J_T.col(2) = [c,d].J_T.col(2)
                            ///     J_T.col(2).[[a,b],[c,d]].[[1],[0]] = J_T.col(2).[a, c]
                            ///     J_T.col(2).[[a,b],[c,d]].[[0],[1]] = J_T.col(2).[b, d]
                            ///     J_T.col(2).[[a,b],[c,d]].J_T.col(2)
                        }
                    }
                }
            }

            /// find out which grid matched best
            max_idx = 0;
            max_score = std::numeric_limits<double>::lowest();
            for(std::size_t i = 0 ; i < 4 ; ++i) {
                if(score[i] > max_score) {
                    max_score = score[i];
                    max_idx = i;
                }
            }

            /// compute the hessian with the result
            HessianType  &hessian_entry = hessian[max_idx];
            GradientType &gradient_entry = gradient[max_idx];

            /// insert positive definite gurantee here
            double off = hessian_entry.maxCoeff() - hessian_entry.minCoeff();
            for(std::size_t i = 0 ; i < 3 ; ++i) {
                hessian_entry(i,i) += off;
            }

            /// solve equeation here
            delta_p = GradientType::Zero();
            delta_p = hessian_entry.fullPivLu().solve(gradient_entry);
            tx  += delta_p(0);
            ty  += delta_p(1);
            phi += delta_p(2);


            /// check for convergence
            if((epsTrans(tx, tx_old) &&
                epsTrans(ty, ty_old) &&
                epsRot(phi, phi_old)))
                break;

            ++iteration;
            if(iteration >= params.max_iterations)
                break;
        }

        return max_score;
    }
private:
    typename GridType::Ptr grid;
};
}
}

#endif // MATCHER2D_HPP
