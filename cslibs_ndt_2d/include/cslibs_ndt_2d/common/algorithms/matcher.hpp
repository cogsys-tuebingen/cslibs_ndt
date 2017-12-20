#ifndef CSLIBS_NDT_2D_MATCHER_HPP
#define CSLIBS_NDT_2D_MATCHER_HPP

#include <cslibs_math_2d/linear/pointcloud.hpp>

#include <cslibs_time/stamped.hpp>
#include <cslibs_gridmaps/static_maps/algorithms/normalize.hpp>
#include <cslibs_gridmaps/static_maps/conversion/convert_probability_gridmap.hpp>
#include <cslibs_gridmaps/utility/delegate.hpp>
#include <nav_msgs/OccupancyGrid.h>

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

    using callback_t            = delegate<void(const nav_msgs::OccupancyGrid::Ptr &)>;

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

    Matcher() = delete;

    Matcher(const Parameters & params = Parameters()) :
        params_(params),
        callback_([](const nav_msgs::OccupancyGrid::Ptr &) {})
    { }

    inline void setCallback(
            const callback_t & cb)
    {
        callback_ = cb;
    }

    inline double match(
            const cslibs_math_2d::Pointcloud2d::Ptr & dst,
            const cslibs_math_2d::Pointcloud2d::Ptr & src,
            transform_t                             & transformation,
            const transform_t                       & prior_transformation = transform_t::identity())
    const
    {
        typename map_t::Ptr dst_map = toMap(dst, transform_t::identity(), params_.resolution_);
        return match(dst_map, src, transformation, prior_transformation);
    }

    inline double match(
            const typename map_t::Ptr               & dst,
            const cslibs_math_2d::Pointcloud2d::Ptr & src,
            transform_t                             & transformation,
            const transform_t                       & prior_transformation = transform_t::identity())
    const
    {
        publish(dst);

        auto wrap_angle = [](double x) {
            x = fmod(x + M_PI, 2.0 * M_PI);
            if (x < 0)
                x += 2.0 * M_PI;
            return x - M_PI;
        };

        std::size_t iterations       = 0;
        std::size_t step_corrections = 0;
        lambda_t    lambda           = params_.lambda_;
        double      max_score        = std::numeric_limits<double>::lowest();

        transform_t correction       = transform_t::identity();
        transform_t best_correction  = transform_t::identity();
        transformation               = prior_transformation * correction;//correction * prior_transformation;

        double score = 0.0;
        Eigen::Vector3d prev_delta;
        for (; iterations < params_.max_iterations_ && step_corrections < params_.max_step_corrections_;
             ++ iterations) {

            gradient_t gradient(gradient_t::Zero());
            hessian_t  hessian (hessian_t::Zero());

            derivatives(dst, src, transformation, hessian, gradient);
            ensurePositiveDefiniteness(hessian, gradient);
/*
            // adapt step size
            if (score > 0.0 && score < max_score) {
                correction     = best_correction;
                lambda        *= params_.alpha_;
                transformation = correction * prior_transformation;
                ++ step_corrections;
            } else {
                if (iterations > 0 && eps(prev_delta))
                    break;

                max_score        = score;
                best_correction  = correction;
                step_corrections = 0;
                lambda          /= params_.alpha_;
            }
//*/
            Eigen::Vector3d delta = hessian.ldlt().solve(-gradient);
            //for (int i = 0 ; i < 3 ; ++ i)
            //    delta(i) *= lambda(i);
            delta = -delta;
            prev_delta = delta;

            // update correction
            Eigen::Vector3d corr = delta + correction.toEigen();
            corr(3) = wrap_angle(corr(3));
            correction.setFrom(corr);

            // update transformation
            transformation = prior_transformation * correction;//correction * prior_transformation;

            // estimate score
            score = scoreCloud(dst, src, transformation);//*
            if (true) {//(score > max_score) {
                best_correction = correction;
                max_score       = score;
            }

            // check epsilon criterion
            if (eps(delta))
                break;//*/
        }

        transformation = prior_transformation * best_correction;//best_correction * prior_transformation;
        return max_score;
    }

    inline double match(
            const typename map_t::Ptr & dst,
            const typename map_t::Ptr & src,
            transform_t               & transformation,
            const transform_t         & prior_transformation = transform_t::identity())
    const
    {
        // TODO: implement
        std::cerr << "not implemented yet" << std::endl;

        return 0.0;
    }

protected:
    virtual typename map_t::Ptr toMap(
            const cslibs_math_2d::Pointcloud2d::Ptr & cloud,
            const pose_t                            & origin,
            const double                            & resolution)
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

    inline bool eps(
            const Eigen::Vector3d & diff)
    const
    {
        return eps(diff(0), diff(1), diff(2));
    }

    inline void derivatives(
            const typename map_t::Ptr               & dst,
            const cslibs_math_2d::Pointcloud2d::Ptr & src,
            const transform_t                       & transformation,
            hessian_t                               & hessian,
            gradient_t                              & gradient)
    const
    {
        const double bundle_resolution = dst->getBundleResolution();
        auto to_bundle_index = [&bundle_resolution] (const point_t & p) -> typename map_t::index_t {
            return {{static_cast<int>(std::floor(p(0) * bundle_resolution)),
                            static_cast<int>(std::floor(p(1) * bundle_resolution))}};
        };

        auto sample = [](const distribution_t * d,
                const point_t        & p,
                point_t              & q) {
            return d ? d->data().sampleNonNormalized(p, q) : 0.0;
        };

        for (const point_t & p : *src) {
            if (!p.isNormal())
                continue;

            const point_t point = transformation * p;
            if (!point.isNormal())
                continue;

            const typename map_t::index_t bi = to_bundle_index(point);
            const distribution_bundle_t * b  = dst->getDistributionBundle(bi);
            if (!b)
                continue;

            int     idx   = -1;
            double  score = 0.0;
            point_t point_mean;
/*/
            for (int i = 0 ; i < 4 ; ++ i) {
                point_t q_tmp;
                const double s = sample(b->at(i), point, q_tmp);
                if (s > score) {
                    idx        = i;
                    score      = s;
                    point_mean = q_tmp;
                }
            }
            if (idx < 0 || score == 0.0)
                continue;
           // std::cout << score << ", " << b->at(idx)->data().sample(point) << std::endl;
/*/
            for (int idx = 0 ; idx < 4 ; ++ idx) {
                point_t point_mean;
                const double score = sample(b->at(idx), point, point_mean);
                if (score == 0.0)
                    continue;//*/

                // first and second order derivatives
                Eigen::Matrix<double, 2, 3> jac(Eigen::Matrix<double, 2, 3>::Identity());
                jac(0, 2) = -point(0) * transformation.sin() - point(1) * transformation.cos();
                jac(1, 2) =  point(0) * transformation.cos() - point(1) * transformation.sin();

                Eigen::Vector2d hes;
                hes(0) = -point(0) * transformation.cos() + point(1) * transformation.sin();
                hes(1) = -point(0) * transformation.sin() - point(1) * transformation.cos();

                // information matrix of the distribution
                const typename distribution_t::distribution_t & d = b->at(idx)->data();
                const Eigen::Matrix2d & information_matrix        = d.getInformationMatrix();

                // use score of q as weighting factor
                Eigen::Vector2d q; q << point_mean(0), point_mean(1);
                const double factor = score / static_cast<double>(src->size());

                // gradient
                gradient_t gradient_tmp(gradient_t::Zero());
                for (int i = 0 ; i < 3 ; ++ i) {
                    gradient_tmp(i) = q.dot(information_matrix * jac.col(i));
                    gradient(i)    -= gradient_tmp(i) * factor;
                }

                // hessian
                for (int i = 0 ; i < 3 ; ++ i) {
                    for (int j = 0  ; j < 3 ; ++ j) {
                        Eigen::Vector2d hij = (i == 2 && j == 2) ? hes : Eigen::Vector2d::Zero();
                        hessian(i, j) -= factor * (-gradient_tmp(i) * gradient_tmp(j) +
                                                   q.dot(information_matrix * hij) +
                                                   jac.col(j).dot(information_matrix * jac.col(i)));
                    }
                }
            }
        }
    }

    inline void ensurePositiveDefiniteness(
            hessian_t        & hessian,
            const gradient_t & gradient)
    const
    {
        Eigen::SelfAdjointEigenSolver<hessian_t> hessian_solver(hessian);
        gradient_t evals = hessian_solver.eigenvalues().real();
        double min_coeff = evals.minCoeff();
        double max_coeff = evals.maxCoeff();

        if(min_coeff < 0) {
            hessian_t evecs = hessian_solver.eigenvectors().real();
            double regularizer = gradient.norm();
            regularizer = regularizer + min_coeff > 0 ? regularizer : 0.001 * max_coeff - min_coeff;

            gradient_t reg;
            reg<<regularizer,regularizer,regularizer;
            evals += reg;

            hessian_t lambda_hessian;
            lambda_hessian = evals.asDiagonal();
            hessian = evecs*lambda_hessian*(evecs.transpose());
        }
    }

    inline double scoreCloud(
            const typename map_t::Ptr               & dst,
            const cslibs_math_2d::Pointcloud2d::Ptr & src,
            const transform_t                       & transformation)
    const
    {
        const double bundle_resolution = dst->getBundleResolution();
        auto to_bundle_index = [&bundle_resolution] (const point_t & p) -> typename map_t::index_t {
            return {{static_cast<int>(std::floor(p(0) * bundle_resolution)),
                            static_cast<int>(std::floor(p(1) * bundle_resolution))}};
        };

        auto sample = [](const distribution_t * d,
                         const point_t        & p) {
            return d ? d->data().sampleNonNormalized(p) : 0.0;
        };
        auto sample_bundle = [&sample] (const distribution_bundle_t * b,
                const point_t               & p)
        {
            return 0.25 * (sample(b->at(0), p) +
                           sample(b->at(1), p) +
                           sample(b->at(2), p) +
                           sample(b->at(3), p));
        };
        auto print = [](const distribution_t * d, const point_t & p) {
            return d ? (std::to_string(d->data().getN()) + ": " + std::to_string(d->data().sample(p))) : "-";
        };

        double score = 0.0;
        for (const point_t & p : *src) {
            if (!p.isNormal())
                continue;

            const point_t point = transformation * p;
            if (!point.isNormal())
                continue;

            const typename map_t::index_t bi = to_bundle_index(point);
            const distribution_bundle_t * b  = dst->getDistributionBundle(bi);
            if (b)
                score += sample_bundle(b, point);
/*
            std::cout << point << " -> " << sample_bundle(b, point) << "; "
                      << print(b->at(0), point) << "; "
                      << print(b->at(1), point) << "; "
                      << print(b->at(2), point) << "; "
                      << print(b->at(3), point) << std::endl;*/
        }

        return score;
    }

    inline void publish(
            const typename map_t::Ptr & dst)
    const
    {
        if (!dst)
            return;

        using static_map_t         = cslibs_gridmaps::static_maps::ProbabilityGridmap;
        using static_map_stamped_t = cslibs_time::Stamped<typename static_map_t::Ptr>;

        const double sampling_resolution_ = 0.025;
        const std::size_t height = static_cast<std::size_t>(dst->getHeight() / sampling_resolution_);
        const std::size_t width  = static_cast<std::size_t>(dst->getWidth()  / sampling_resolution_);
        std::cout << dst->getResolution() << ", " << height << ", " << width << std::endl;

        cslibs_math_2d::Transform2d origin = dst->getOrigin();
        static_map_stamped_t map_tmp;
        map_tmp.data().reset(new static_map_t(origin,
                                              sampling_resolution_,
                                              height,
                                              width));
std::cout << origin << std::endl;
        const double bundle_resolution = dst->getBundleResolution();
        const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution_);
        typename map_t::index_t min_distribution_index = dst->getMinDistributionIndex();
        typename map_t::index_t max_distribution_index = dst->getMaxDistributionIndex();

        auto sample = [](const typename map_t::distribution_t *d,
                         const cslibs_math_2d::Point2d &p) {
            return d ? d->data().sampleNonNormalized(p)/*/(d->data().getN() > 0 ? 1.0 : 0.0)*/ : 0.0;
        };
        auto sample_bundle = [&sample] (const typename map_t::distribution_bundle_t* b,
                                        const cslibs_math_2d::Point2d &p)
        {
            return 0.25 * (sample(b->at(0), p) +
                           sample(b->at(1), p) +
                           sample(b->at(2), p) +
                           sample(b->at(3), p));
        };

        for(int i = min_distribution_index[1] ; i <= max_distribution_index[1] ; ++i)
            for(int j = min_distribution_index[0] ; j <= max_distribution_index[0] ; ++j)
                dst->getDistributionBundle({{j,i}});

        min_distribution_index = dst->getMinDistributionIndex();
        max_distribution_index = dst->getMaxDistributionIndex();

        for(int i = min_distribution_index[1] ; i <= max_distribution_index[1] ; ++i) {
            for(int j = min_distribution_index[0] ; j <= max_distribution_index[0] ; ++j) {
                const typename map_t::distribution_bundle_t* bundle = dst->getDistributionBundle({{j,i}});
                if(bundle) {
                    const int cx = (j - min_distribution_index[0]) * static_cast<int>(chunk_step);
                    const int cy = (i - min_distribution_index[1]) * static_cast<int>(chunk_step);
                    for(int k = 0 ; k < chunk_step ; ++k) {
                        for(int l = 0 ; l < chunk_step ; ++l) {
                            const cslibs_math_2d::Point2d p(j * bundle_resolution + l * sampling_resolution_,
                                                            i * bundle_resolution + k * sampling_resolution_);
                            map_tmp.data()->at(static_cast<std::size_t>(cx + l),
                                               static_cast<std::size_t>(cy + k)) = sample_bundle(bundle, p);
                        }
                    }
                }
            }
        }
        cslibs_gridmaps::static_maps::algorithms::normalize<double>(*map_tmp.data());

        nav_msgs::OccupancyGrid::Ptr map;
        cslibs_gridmaps::static_maps::conversion::from(map_tmp, map);
        map->header.stamp    = ros::Time::now();
        map->header.frame_id = "laser";
        callback_(map);
    }

private:
    Parameters params_;
    callback_t callback_;

};
}
}
}

#endif // CSLIBS_NDT_2D_MATCHER_HPP
