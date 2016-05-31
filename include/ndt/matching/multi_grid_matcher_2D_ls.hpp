#ifndef MATCHER2DLS_HPP
#define MATCHER2DLS_HPP

#include <ndt/grid/multi_grid.hpp>
#include <ndt/matching/matcher.hpp>
#include <ndt/data/pointcloud.hpp>
#include <ndt/math/angle.hpp>
#include <ndt/math/hausdorff.hpp>

#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace ndt {
namespace matching {
class MultiGridMatcher2DLS : public Matcher<2> {
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

    struct LineSearch {
        double min;
        double max;
        double value;
        bool   lower;

        LineSearch(const double _min,
                   const double _max) :
            min(_min),
            max(_max),
            value((max + min) / 2.0),
            lower(false)
        {
        }

        void reset(const double _min,
                   const double _max)
        {
            min = _min;
            max = _max;
            value = ((max + min)  / 2.0);
            lower = false;
        }

        double next(const bool success)
        {
            if(success) {
                if(lower) {
                    max = value;
                } else {
                    min = value;
                }
                value = (max + min) / 2.0;
            } else {
                lower = !lower;
            }

            if(lower) {
                return (value + min) / 2.0;
            } else {
                return (max + value) / 2.0;
            }
        }



    };


    MultiGridMatcher2DLS(const Parameters &_params = Parameters()) :
        BaseClass(_params),
        rotation(0.0),
        line_search(0.001, _params.alpha)
    {
    }

    virtual ~MultiGridMatcher2DLS()
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
        tx        = 0.0;
        ty        = 0.0;
        phi       = 0.0;
        prev_tx  = tx;
        prev_ty  = ty;
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
            if(iteration > 0 &&
                    epsTrans(tx, prev_tx) &&
                    epsTrans(ty, prev_ty) &&
                    epsRot(phi, prev_phi)) {
                break;
            }

            prev_tx   = tx;
            prev_ty   = ty;
            prev_phi  = phi;

            rotation        = RotationType(phi);
            translation     = TranslationType(tx, ty);
            transformation  = translation * rotation * _prior_transformation;

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

                        double s = distribution.sampleNonNormalized(p, q);
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

            if(iteration >= params.max_iterations) {
                break;
            }

            std::size_t current_max_idx = 0;
            double current_max_score = std::numeric_limits<double>::lowest();
            for(std::size_t i = 0 ; i < 4 ; ++i) {
                if(score[i] > current_max_score) {
                    current_max_score = score[i];
                    current_max_idx = i;
                }
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

            /// here we apply the line search
            max_score = lineSearch(_src, max_score, delta_p, _prior_transformation);
            ++iteration;
        }

        rotation        = RotationType(phi);
        translation     = TranslationType(tx, ty);
        transformation  = translation * rotation * _prior_transformation;
        _transformation = transformation;

        /// todo exchange through kdtree nn ...
        std::size_t accepted = 0;
        std::size_t size_valid = 0;
        for(std::size_t i = 0 ; i < _dst.size ; ++i) {
            if(_src.mask[i] == PointCloudType::VALID) {
                PointType p = transformation * _src.points[i];
                double h = ndt::math::hausdorff<2>(p, _dst);
                if(h < 0.5)
                    ++accepted;
                ++size_valid;
            }
        }

        return accepted / (double) size_valid;

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

protected:
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
    LambdaType                  lambda;
    LineSearch                  line_search;
    std::size_t                 step_corrections;

    double lineSearch(const PointCloudType   &_src,
                      const double         &_score,
                      const GradientType   &_delta_p,
                      const TransformType  &_prior_transformation)
    {

        line_search.reset(0.001, params.alpha);

        double ls_tx;
        double ls_ty;
        double ls_phi;
        RotationType    ls_rotation(0.0);
        TranslationType ls_translation;
        TransformType   ls_transformation;
        double ls_tx_max = 0;
        double ls_ty_max = 0;
        double ls_phi_max = 0;
        double max_score = _score;
        double ls_scale  = line_search.value;
        std::array<double, 4>  ls_score;

        bool converged = false;
        std::size_t iteration = 0;
        while(!converged) {
            if(iteration > 100)
                break;

            ls_score.fill(0.0);
            ls_tx  = tx  + _delta_p(0) * ls_scale * lambda(0);
            ls_ty  = ty  + _delta_p(1) * ls_scale * lambda(1);
            ls_phi = phi + _delta_p(2) * ls_scale * lambda(2);

            ls_rotation    = RotationType(ls_phi);
            ls_translation = TranslationType(ls_tx, ls_ty);
            ls_transformation = translation * rotation * _prior_transformation;

            /// calculate the hessian and the gradient for each grid
            for(std::size_t i = 0 ; i < _src.size ; ++i) {
                if(_src.mask[i] == PointCloudType::VALID) {
                    PointType p = ls_transformation * _src.points[i];

                    DistributionsType distributions;
                    grid->get(p, distributions);
                    for(std::size_t j = 0 ; j < distributions.size(); ++j) {
                        if(distributions[j] == nullptr)
                            continue;
                        DistributionType &distribution = *distributions[j];
                        if(distribution.getN() < 3)
                            continue;

                        double s = distribution.sampleNonNormalized(p);
                        if(s > 1e-3) {
                            ls_score[j] += s;
                        }
                    }
                }
            }

            /// find out which grid matched best
            double current_max_score = std::numeric_limits<double>::lowest();
            for(std::size_t i = 0 ; i < 4 ; ++i) {
                if(score[i] > current_max_score) {
                    current_max_score = ls_score[i];
                }
            }

            std::cout << max_score << " " << current_max_score << std::endl;
            std::cout << ls_scale << std::endl;


            if(eps(max_score, current_max_score, 1e-2))
                break;

            if(current_max_score >= max_score) {
                ls_scale = line_search.next(true);
                ls_tx_max   = ls_tx;
                ls_ty_max   = ls_ty;
                ls_phi_max = ls_phi;
                max_score = current_max_score;
            } else {
                ls_scale = line_search.next(false);
            }
            ++iteration;
        }

        tx  = ls_tx_max;
        ty  = ls_ty_max;
        phi = math::wrapAngle(ls_phi_max);
        return max_score;

    }

};
}
}

#endif // MATCHER2DLS_HPP
