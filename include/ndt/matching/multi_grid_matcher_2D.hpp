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
    typedef Eigen::Matrix<double,2,3>                  JacobiType;
    typedef Eigen::Translation2d                       TranslationType;
    typedef Eigen::Rotation2Dd                         RotationType;
    typedef Eigen::Vector3d                            GradientType;


    MultiGridMatcher2D(const Parameters &_params = Parameters()) :
        BaseClass(_params),
        rotation(0.0)
    {
    }

    virtual ~MultiGridMatcher2D()
    {
    }

    inline double match(const PointCloudType &_dst,
                        const PointCloudType &_src,
                        TransformType        &_transformation,
                        const TransformType  &_prior_transformation = TransformType::Identity()) override
    {
        /// build the ndt grid for the src cloud
        PointType range = _dst.range();
        typename GridType::SizeType size;
        for(std::size_t i = 0 ; i < 2 ; ++i) {
            if(range(i) <= 0.0)
                throw std::range_error("Point cloud boundaries are not set properly!");
            size[i] = floor(range(i) / params.resolution[i] + 0.5);
        }

        grid.reset(new GridType(size, params.resolution, _dst.min));
        grid->add(_dst);

        /// reset all the members
        tx        = _prior_transformation.translation()(0);
        ty        = _prior_transformation.translation()(1);
        phi       = atan2(_prior_transformation.rotation()(1,0),
                          _prior_transformation.rotation()(0,0));
        prev_tx = tx;
        prev_ty = ty;
        prev_phi = phi;
        /// gradient and stuff
        /// need 4 fields for that
        /// only solve the maximal score
        max_score             = std::numeric_limits<double>::lowest();
        iteration             = 0;
        lambda                = params.lambda;
        step_corrections      = 0;



        bool converged = false;
        while(!converged) {
            rotation        = RotationType(phi);
            translation     = TranslationType(tx, ty);
            transformation  = translation * rotation;

            gradient.fill(GradientType::Zero());
            hessian.fill(HessianType::Zero());
            score.fill(0.0);
            double  sin_phi;
            double  cos_phi;
            sincos(phi, &sin_phi, &cos_phi);

            /// calculate the hessian and the gradient for each grid
            for(std::size_t i = 0 ; i < _src.size ; ++i) {
                if(_src.mask[i] == PointCloudType::VALID) {
                    PointType p = transformation * _src.points[i];

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
                        PointType            jac = PointType::Zero();
                        PointType            hes = PointType::Zero();
                        GradientType        &gradient_entry = gradient[j];
                        HessianType         &hessian_entry  = hessian[j];
                        double              &score_entry = score[j];

                        double s = distribution.sampleNonNoramlized(p, q);
                        distribution.getMean(mean);
                        distribution.getInverseCovariance(inverse_covariance);

                        /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                        /// if s is 0.0 we do no need to count the sample in
                        if(s > 1e-3) {
                            score_entry += s;
                            const PointType q_inverse_covariance = q.transpose() * inverse_covariance;

                            jac(0) = -p(0) * sin_phi - p(1) * cos_phi;
                            jac(1) =  p(0) * cos_phi - p(1) * sin_phi;
                            hes(0) = -p(0) * cos_phi + p(1) * sin_phi;
                            hes(1) = -p(0) * sin_phi - p(1) * cos_phi;

                            //                            /// gradient computation
                            double g_dot = q_inverse_covariance.transpose() * jac;
                            gradient_entry(0) += s * q_inverse_covariance(0);
                            gradient_entry(1) += s * q_inverse_covariance(1);
                            gradient_entry(2) += s * g_dot;

                            hessian_entry(0,0)+= -s * (
                                        (-q_inverse_covariance(0) * -q_inverse_covariance(0))   /// (1)
                                        +(-inverse_covariance(0,0)));                           /// (3)
                            hessian_entry(1,0)+= -s * (
                                        (-q_inverse_covariance(1) * -q_inverse_covariance(0))   /// (1)
                                        +(-inverse_covariance(1,0)));                           /// (3)
                            hessian_entry(2,0)+= -s * (
                                        (-g_dot * -q_inverse_covariance(0))                     /// (1)
                                        +(-inverse_covariance.row(0).dot(jac)));                /// (3)
                            hessian_entry(0,1)+= -s * (
                                        (-q_inverse_covariance(0) * -q_inverse_covariance(1))   /// (1)
                                        +(-inverse_covariance(0,1)));                           /// (3)
                            hessian_entry(1,1)+= -s * (
                                        (-q_inverse_covariance(1) * -q_inverse_covariance(1))   /// (1)
                                        +(-inverse_covariance(1,1)));                           /// (3)
                            hessian_entry(2,1)+= -s * (
                                        (-g_dot * -q_inverse_covariance(1))                     /// (1)
                                        +(-inverse_covariance.row(1).dot(jac)));                /// (3)
                            hessian_entry(0,2)+= -s * (
                                        (-q_inverse_covariance(0) * -g_dot)                     /// (1)
                                       +(-jac.transpose() * inverse_covariance.col(0)));        /// (3)
                            hessian_entry(1,2)+= -s * (
                                        (-q_inverse_covariance(1) * -g_dot )                    /// (1)
                                       +(-jac.transpose() * inverse_covariance.col(1)));        /// (3)
                            hessian_entry(2,2)+= -s * (
                                        (-g_dot * -g_dot)                                       /// (1)
                                        +(-q_inverse_covariance.transpose() * hes)              /// (2)
                                        +(-jac.transpose() * inverse_covariance * jac));        /// (3)


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
            std::size_t current_max_idx = 0;
            double current_max_score = std::numeric_limits<double>::lowest();
            for(std::size_t i = 0 ; i < 4 ; ++i) {
                if(score[i] > current_max_score) {
                    current_max_score = score[i];
                    current_max_idx = i;
                }
            }
            /// now we have to check wether the score increased or not
            /// if not, we have to adjust the step size
            std::cout << "curr and max : " << current_max_score << " " << max_score << std::endl;
            if(current_max_score < max_score) {
                tx      = prev_tx;
                ty      = prev_ty;
                phi     = prev_phi;
                lambda *= params.alpha;
                rotation        = RotationType(phi);
                translation     = TranslationType(tx, ty);
                transformation  = translation * rotation;
                ++step_corrections;
            } else {
                if(iteration > 0 &&
                        epsTrans(tx, prev_tx) &&
                        epsTrans(ty, prev_ty) &&
                        epsRot(phi, prev_phi)) {
                    std::cout << "eps converged" << std::endl;
                    break;
                }

                max_score = current_max_score;
                prev_tx   = tx;
                prev_ty   = ty;
                prev_phi  = phi;
                step_corrections = 0;
            }

            if(step_corrections >= params.max_step_corrections) {
                std::cout << "too many step corrections : " << step_corrections << " vs. " << params.max_step_corrections << std::endl;
                break;
            }

            if(iteration >= params.max_iterations) {
                std::cout << "reached max iterations" <<  std::endl;
                break;
            }


            /// compute the hessian with the result
            HessianType  &hessian_entry  = hessian[current_max_idx];
            GradientType &gradient_entry = gradient[current_max_idx];
            /// positive definiteness
            Eigen::EigenSolver<HessianType> hessian_solver(hessian_entry, false);
            double min_eigen_value = hessian_solver.eigenvalues().real().minCoeff();
            if(min_eigen_value < 0) {
                double lambda = 1.1 * min_eigen_value - 1;
                hessian_entry += Eigen::Vector3d (-lambda, -lambda, -lambda).asDiagonal();
            }

            /// solve equeation here
            delta_p = GradientType::Zero();
            delta_p = (-hessian_entry).fullPivLu().solve(gradient_entry);
            tx  += delta_p(0) * lambda;
            ty  += delta_p(1) * lambda;
            phi += delta_p(2) * lambda;

            ++iteration;
        }

        rotation        = RotationType(phi);
        translation     = TranslationType(tx, ty);
        transformation  = translation * rotation;
        _transformation = transformation;

        return max_score;
    }

    void printDebugInfo()
    {
        std::cout << "----------------------------------" << std::endl;
        std::cout << "score : " << max_score << std::endl;
        std::cout << "translation : " << std::endl;
        std::cout << transformation.translation() << std::endl;
        std::cout << "rotation : " << std::endl;
        std::cout << transformation.rotation() << std::endl;
        std::cout << "iterations : " << iteration << std::endl;
        std::cout << "lambda : " << lambda << std::endl;
        std::cout << "step corrections : " << step_corrections << std::endl;
        std::cout << "----------------------------------" << std::endl;

    }

private:
    typename GridType::Ptr grid;

    double tx;
    double ty;
    double phi;
    double prev_tx;
    double prev_ty;
    double prev_phi;
    double max_score;

    TransformType               transformation;
    TranslationType             translation;
    RotationType                rotation;

    std::array<double, 4>       score;
    std::array<GradientType, 4> gradient;
    std::array<HessianType, 4>  hessian;
    GradientType                delta_p;

    std::size_t                 iteration;
    double                      lambda;
    std::size_t                 step_corrections;


};
}
}

#endif // MATCHER2D_HPP
