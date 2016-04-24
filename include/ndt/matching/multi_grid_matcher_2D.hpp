
#ifndef MATCHER2D_HPP
#define MATCHER2D_HPP

#include <ndt/grid/multi_grid.hpp>
#include <ndt/matching/matcher.hpp>
#include <ndt/data/pointcloud.hpp>

#include <fstream>
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace ndt {
namespace matching {
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class MultiGridMatcher2D : public NDTMatcher<2> {
public:
    typedef Eigen::Matrix<double,3,3>  HessianType;
    typedef Eigen::Translation2d       TranslationType;
    typedef Eigen::Rotation2Dd         RotationType;

    MultiGridMatcher2D(const ResolutionType &resolution) :
        BaseClass(resolution),
        out("/tmp/out.txt")
    {
    }

    virtual ~MultiGridMatcher2D()
    {
        out.close();
    }

private:
    std::ofstream out;

    double doMatch(const PointCloudType &_dst,
                   TransformType        &_transformation,
                   const std::size_t     _max_iterations = 100,
                   const double          _eps = 1e-3,
                   const double          _eps_rad = 1e-6) override
    {
        _transformation.setIdentity();
        /// todo:: initialize parameter estimate double phi = 0.0;
        double tx = 0.0;
        double ty = 0.0;
        double phi = 0.0;
        TranslationType trans;
        RotationType    rotation(0.0);

        /// variables needed for sampling a point
        PointType            mean;
        CovarianceMatrixType inverse_covariance;
        PointType            q;
        double               s;
        double               g_dot;
        PointType            q_inverse_covariance;
        double               sin_phi;
        double               cos_phi;
        PointType            jac;
        PointType            hes;


        /// gradient and stuff
        GradientType gradient;
        GradientType delta_p;
        HessianType  hessian;
        double       score;
        double       tx_old;
        double       ty_old;
        double       phi_old;
        std::size_t  iteration = 0;

        bool converged = false;
        while(!converged) {
            rotation = RotationType(phi);
            trans    = TranslationType(tx, ty);
            _transformation = trans * rotation;

            gradient = GradientType::Zero();
            hessian  = HessianType::Zero();
            score    = 0.0;
            tx_old   = tx;
            ty_old   = ty;
            phi_old  = phi;
            sincos(phi, &sin_phi, &cos_phi);

            for(std::size_t i = 0 ; i < _dst.size ; ++i) {
                if(_dst.mask[i] == PointCloudType::VALID) {
                    PointType p = _transformation * _dst.points[i];

                    DistributionsType distributions;
                    grid->get(p, distributions);
                    for(std::size_t j = 0 ; j < distributions.size(); ++j) {
                        DistributionType &distribution = *distributions[j];
                        if(distribution.getN() < 3)
                            continue;

                        s = distribution.evaluateNonNoramlized(p, q);
                        distribution.getMean(mean);
                        distribution.getInverseCovariance(inverse_covariance);

                        /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                        /// if s is 0.0 we do no need to count the sample in
                        if(s > 0.0) {
                            score += s;
                            q_inverse_covariance = q.transpose() * inverse_covariance;

                            jac(0) = -q(0) * sin_phi - q(1) * cos_phi;
                            jac(1) =  q(0) * cos_phi - q(1) * sin_phi;
                            hes(0) = -q(0) * cos_phi + q(1) * sin_phi;
                            hes(1) = -q(0) * sin_phi - q(1) * cos_phi;

                            /// gradient computation
                            g_dot = q_inverse_covariance.dot(jac);
                            gradient(0) -= s * q_inverse_covariance(0);
                            gradient(1) -= s * q_inverse_covariance(1);
                            gradient(2) -= s * g_dot;
                            /// hessian computation
                            hessian(0,0)+=  s
                                    * (-q_inverse_covariance(0) * q_inverse_covariance(0)   /// (1)
                                       +  inverse_covariance(0,0));                            /// (3)
                            hessian(1,0)+=  s
                                    * (-q_inverse_covariance(1) * q_inverse_covariance(0)   /// (1)
                                       +  inverse_covariance(1,0));                            /// (3)
                            hessian(2,0)+=  s
                                    * (-g_dot * inverse_covariance(0)                       /// (1)
                                       +  inverse_covariance.row(0) * jac);                    /// (3)
                            hessian(0,1)+=  s
                                    * (-q_inverse_covariance(0) * q_inverse_covariance(1)   /// (1)
                                       +  inverse_covariance(0,1));                            /// (3)
                            hessian(1,1)+=  s
                                    * (-q_inverse_covariance(1) * q_inverse_covariance(1)   /// (1)
                                       +  inverse_covariance(1,1));                            /// (3)
                            hessian(2,1)+=  s
                                    * (-g_dot * inverse_covariance(1)                       /// (1)
                                       +  inverse_covariance.row(1) * jac);                    /// (3)
                            hessian(0,2)+=  s
                                    * (-q_inverse_covariance(0) * g_dot                     /// (1)
                                       +  jac.transpose() * inverse_covariance.col(0));        /// (3)
                            hessian(1,2)+=  s
                                    * (-q_inverse_covariance(1) * g_dot                     /// (1)
                                       +  jac.transpose() * inverse_covariance.col(1));        /// (3)
                            hessian(2,2)+=  s
                                    * (-g_dot * g_dot                                       /// (1)
                                       +  q_inverse_covariance.dot(hes)                        /// (2)
                                       +  jac.transpose() * inverse_covariance * jac);         /// (3)

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
            /// insert positive definite gurantee here
            double max = std::numeric_limits<double>::min();
            double min = std::numeric_limits<double>::max();
            for(std::size_t i = 0 ; i < 3 ; ++i) {
                if(hessian(i,i) > max)
                    max = hessian(i,i);
                if(hessian(i,i) < min)
                    min = hessian(i,i);
            }
            double off = max - min;
            for(std::size_t i = 0 ; i < 3 ; ++i) {
                hessian(i,i) += off;
            }

            /// solve equeation here
            delta_p = GradientType::Zero();
            delta_p = hessian.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);
            tx  += delta_p(0);
            ty  += delta_p(1);
            phi += delta_p(2);

            //            out << tx << " " << ty << " " << phi << std::endl;

            /// check for convergence
            if((eps(tx, tx_old, _eps) &&
                eps(ty, ty_old, _eps) &&
                eps(phi, phi_old, _eps_rad)) ||
                    iteration > _max_iterations)
                break;
            ++iteration;
        }

//        std::cout << "iterations " << iteration << " " << score << std::endl;

        return score;

    }

    inline bool eps(const double a,
                    const double b,
                    const double epsilon = 1e-3)
    {
        return fabs(a - b) < epsilon;
    }
};
}
}

#endif // MATCHER2D_HPP
