#ifndef CSLIBS_NDT_2D_MATCHER_HPP
#define CSLIBS_NDT_2D_MATCHER_HPP

#include <cslibs_math_2d/linear/pointcloud.hpp>
#include <cslibs_math/common/floor.hpp>

namespace cslibs_ndt_2d {
namespace common {
namespace algorithms {
template <typename map_t>
class Matcher {
public:
    using transform_t           = typename map_t::transform_t;
    using point_t               = typename map_t::point_t;
    using pose_t                = typename map_t::pose_t;
    using distribution_t        = typename map_t::distribution_t;
    using distribution_bundle_t = typename map_t::distribution_bundle_t;

    using lambda_t              = Eigen::Vector3d;
    using gradient_t            = Eigen::Vector3d;
    using hessian_t             = Eigen::Matrix3d;

    struct Parameters {
        Parameters() :
            resolution_(1.0),
            eps_rot_(1e-3),
            eps_trans_(1e-3),
            max_iterations_(100),
            max_step_corrections_(10),
            lambda_(lambda_t::Constant(0.1)),
            alpha_(2.0)
        { }

        Parameters(
                const double      & resolution,
                const double      & eps_rot,
                const double      & eps_trans,
                const std::size_t & max_iterations,
                const std::size_t & max_step_corrections,
                const lambda_t    & lambda,
                const double      & alpha) :
            resolution_(resolution),
            eps_rot_(eps_rot),
            eps_trans_(eps_trans),
            max_iterations_(max_iterations),
            max_step_corrections_(max_step_corrections),
            lambda_(lambda),
            alpha_(alpha)
        { }

        double      resolution_;
        double      eps_rot_;
        double      eps_trans_;
        std::size_t max_iterations_;
        std::size_t max_step_corrections_;
        lambda_t    lambda_;
        double      alpha_;
    };

    Matcher(const Parameters & params = Parameters()) :
        params_(params)
    { }

    inline double match(
            const cslibs_math_2d::Pointcloud2d & dst,
            const cslibs_math_2d::Pointcloud2d & src,
            transform_t                        & transformation,
            const transform_t                  & prior_transformation = transform_t::identity())
    const
    {
        typename map_t::Ptr dst_map = toMap(dst, transform_t::identity(), params_.resolution_);
        return match(*dst_map, src, transformation, prior_transformation);
    }

    inline double match(
            const map_t                        & dst,
            const cslibs_math_2d::Pointcloud2d & src,
            transform_t                        & transformation,
            const transform_t                  & prior_transformation = transform_t::identity())
    const
    {
        const double bundle_resolution = dst.getBundleResolution();
        auto to_bundle_index = [&bundle_resolution] (const point_t & p) -> typename map_t::index_t {
            return {{static_cast<int>(cslibs_math::common::floor(p(0) * bundle_resolution)),
                            static_cast<int>(cslibs_math::common::floor(p(1) * bundle_resolution))}};
        };

        std::size_t iterations          = 0;
        std::size_t step_corrections    = 0;
        lambda_t    lambda              = params_.lambda_;
        double      max_score           = std::numeric_limits<double>::lowest();

        transform_t correction          = transform_t::identity();
        transform_t previous_correction = transform_t::identity();
        transformation                  = correction * prior_transformation;

        for (; iterations < params_.max_iterations_ && step_corrections < params_.max_step_corrections_;
             ++ iterations) {

            std::array<gradient_t, 4> gradients;
            std::array<hessian_t, 4>  hessians;
            std::array<double, 4>     scores;

            gradients.fill(gradient_t::Zero());
            hessians.fill (hessian_t::Zero());
            scores.fill   (0.0);

            // calculate hessian and gradient
            for (const point_t & p : src) {
                if (p.isNormal()) {

                    // transform point with current transformation
                    const point_t point = transformation * p;
                    if (point.isNormal()) {

                        // get distribution bundle
                        const typename map_t::index_t bi = to_bundle_index(point);
                        const distribution_bundle_t * b  = dst.getDistributionBundle(bi);
                        if (b) {

                            // iterate over all distributions of the bundle
                            for (std::size_t i = 0 ; i < 4 ; ++ i) {
                                if (!b->at(i))
                                    continue;

                                // the i-th distribution
                                typename distribution_t::distribution_t d = b->at(i)->getHandle()->data();
                                if (d.getN() < 3)
                                    continue;

                                // estimate the score of point p in this distribution
                                const double s = d.sampleNonNormalized(point);
                                if (s <= 1e-3)
                                    continue;
                                scores[i] += s;

                                // mean and information matrix of distribution d
                                const Eigen::Vector2d & mean               = d.getMean();
                                const Eigen::Matrix2d & information_matrix = d.getInformationMatrix();

                                // point q and q^t * inf
                                Eigen::Vector2d point_eigen;
                                point_eigen(0) = point(0);
                                point_eigen(1) = point(1);
                                const Eigen::Vector2d             q      = point_eigen - mean;
                                const Eigen::Matrix<double, 1, 2> qt_inf = q.transpose() * information_matrix;

                                // first and second order derivatives of point
                                Eigen::Vector2d jac, hes;
                                jac(0) = -q(0) * transformation.sin() - q(1) * transformation.cos();
                                jac(1) =  q(0) * transformation.cos() - q(1) * transformation.sin();
                                hes(0) = -q(0) * transformation.cos() + q(1) * transformation.sin();
                                hes(1) = -q(0) * transformation.sin() - q(1) * transformation.cos();
                                const Eigen::Matrix<double, 1, 2> jact = jac.transpose();

                                // gradient
                                const double g_dot = qt_inf.dot(jac);
                                gradients[i](0) += s * qt_inf(0);
                                gradients[i](1) += s * qt_inf(1);
                                gradients[i](2) += s * g_dot;

                                // hessian
                                hessians[i](0, 0) -= s * ((qt_inf(0) * qt_inf(0)) - information_matrix(0, 0));
                                hessians[i](1, 0) -= s * ((qt_inf(1) * qt_inf(0)) - information_matrix(1, 0));
                                hessians[i](2, 0) -= s * ((g_dot     * qt_inf(0)) - information_matrix.row(0).dot(jac));
                                hessians[i](0, 1) -= s * ((qt_inf(0) * qt_inf(1)) - information_matrix(0, 1));
                                hessians[i](1, 1) -= s * ((qt_inf(1) * qt_inf(1)) - information_matrix(1, 1));
                                hessians[i](2, 1) -= s * ((g_dot     * qt_inf(1)) - information_matrix.row(1).dot(jac));
                                hessians[i](0, 2) -= s * ((qt_inf(0) * g_dot    ) - (jact.dot(information_matrix.col(0))));
                                hessians[i](1, 2) -= s * ((qt_inf(1) * g_dot    ) - (jact.dot(information_matrix.col(1))));
                                hessians[i](2, 2) -= s * ((g_dot     * g_dot    ) - qt_inf.dot(hes) - (jact * information_matrix * jac));
                            }
                        }
                    }
                }
            }

            // find out which grid matched best
            std::size_t current_max_idx   = 0;
            double      current_max_score = std::numeric_limits<double>::lowest();
            for (std::size_t i = 0 ; i < 4 ; ++ i) {
                if (scores[i] > current_max_score) {
                    current_max_score = scores[i];
                    current_max_idx   = i;
                }
            }

            // adapt step size
            if (current_max_score < max_score) {

                // if score has not improved, reset correction to previous correction and increase step size
                correction     = previous_correction;
                transformation = correction * prior_transformation;
                lambda        *= params_.alpha_;
                ++ step_corrections;

            } else {
                // break, if epsilon criterion was already met
                const Eigen::Vector3d diff = previous_correction.toEigen() - correction.toEigen();
                if (iterations > 0 && eps(diff(0), diff(1), diff(2)))
                    break;

                // if score has improved, decrease step size
                if (current_max_score > max_score)
                    lambda /= params_.alpha_;

                max_score           = current_max_score;
                previous_correction = correction;
                step_corrections    = 0;
            }

            // use hessian and gradient of best grid cells
            gradient_t & gradient = gradients[current_max_idx];
            hessian_t  & hessian  = hessians[current_max_idx];

            // ensure positive definiteness of hessian
            Eigen::EigenSolver<hessian_t> hessian_solver(hessian, false);
            const double min_eigen_value = hessian_solver.eigenvalues().real().minCoeff();
            if (min_eigen_value < 0) {
                double lambda_hessian = 1.1 * std::fabs(min_eigen_value) + 1.0;
                hessian += Eigen::Vector3d(lambda_hessian, lambda_hessian, lambda_hessian).asDiagonal();
            }

            // solve equation
            const Eigen::Vector3d delta = lambda.array() * hessian.fullPivLu().solve(gradient).array();

            // update correction
            correction.setFrom(delta + correction.toEigen());

            // update transformation
            transformation = correction * prior_transformation;

            // check epsilon criterion
            if (eps(delta(0), delta(1), delta(2)))
                break;
        }

        return max_score;
    }

    inline double match(
            const map_t       & dst,
            const map_t       & src,
            transform_t       & transformation,
            const transform_t & prior_transformation = transform_t::identity())
    const
    {
        // TODO: implement
        std::cerr << "not implemented yet" << std::endl;

        return 0.0;
    }

protected:
    virtual typename map_t::Ptr toMap(
            const cslibs_math_2d::Pointcloud2d & cloud,
            const pose_t                       & origin,
            const double                       & resolution)
    const = 0;

    inline bool epsTrans(
            const double & diff)
    const
    {
        return std::fabs(diff) < params_.eps_trans_;
    }

    inline bool epsRot(
            const double & diff)
    const
    {
        return std::fabs(diff) < params_.eps_rot_;
    }

    inline bool eps(
            const double & diff_trans_x,
            const double & diff_trans_y,
            const double & diff_rot)
    const
    {
        return epsTrans(diff_trans_x) && epsTrans(diff_trans_y) && epsRot(diff_rot);
    }

private:
    Parameters params_;
};
}
}
}

#endif // CSLIBS_NDT_2D_MATCHER_HPP
