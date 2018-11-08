#pragma once

#include <cslibs_ndt_3d/matching/gridmap_match_traits.hpp>
#include <cslibs_ndt/matching/match_traits.hpp>
#include <cslibs_ndt/matching/parameter.hpp>
#include <cslibs_ndt/matching/result.hpp>

namespace cslibs_ndt {
namespace matching {

template<typename iterator_t, typename ndt_t, typename traits_t = MatchTraits<ndt_t>>
auto match(const iterator_t& points_begin,
           const iterator_t& points_end,
           const ndt_t& map,
           const Parameter& param,
           const typename ndt_t::transform_t& initial_transform)
-> Result<typename ndt_t::transform_t>
{
    static constexpr int DIMS = traits_t::LINEAR_DIMS + traits_t::ANGULAR_DIMS;
    using point_t             = typename ndt_t::point_t;
    using transform_t         = typename ndt_t::transform_t;
    using result_t            = Result<transform_t>;

    using JacobianCompute = typename traits_t::Jacobian;
    using HessianCompute  = typename traits_t::Hessian;

    using linear_t      = Eigen::Matrix<double, traits_t::LINEAR_DIMS, 1>;
    using angular_t     = Eigen::Matrix<double, traits_t::ANGULAR_DIMS, 1>;
    using gradient_t    = Eigen::Matrix<double, DIMS, 1>;
    using hessian_t     = Eigen::Matrix<double, DIMS, DIMS>;

    // todo: pre transform points, should be externalized or made completely optional...
    std::vector<point_t> points_prime;
    points_prime.reserve(std::distance(points_begin, points_end));
    std::transform(points_begin, points_end, std::back_inserter(points_prime),
                   [&](const point_t& point) { return initial_transform * point; });

    // initialize result
    double max_score        = std::numeric_limits<double>::lowest();
    std::size_t iteration   = 0;

    linear_t  linear    = linear_t::Zero();
    angular_t angular   = angular_t::Zero();

    // initialize state
    linear_t  linear_old    = linear_t::Zero();
    angular_t angular_old   = angular_t::Zero();
    linear_t  linear_delta  = linear_t::Constant(std::numeric_limits<double>::max());
    angular_t angular_delta = angular_t::Constant(std::numeric_limits<double>::max());

    double lambda = 1.0;
    std::size_t step_adjustments = 0;

    // termination criteria
    const auto test_eps = [&]()
    {
        return (linear_delta.array().abs() < param.translationEpsilon()).all()
                && (angular_delta.array().abs() < param.rotationEpsilon()).all();
    };

    const auto test_readjustments = [&]()
    {
        return step_adjustments > 0 && step_adjustments > param.maxStepReadjustments();
    };

    // termination
    const auto terminate = [&](Termination reason)
    {
        return result_t{
            max_score,
                    iteration,
                    traits_t::makeTransform(linear, angular) * initial_transform,
                    reason };
    };

    // iterations
    for (iteration = 0; iteration < param.maxIterations(); ++iteration)
    {
        if (test_readjustments())
            return terminate(Termination::MAX_STEP_READJUSTMENTS);

        const auto t = traits_t::makeTransform(linear, angular);

        JacobianCompute J;
        JacobianCompute::get(angular, J);
        HessianCompute H;
        HessianCompute::get(angular, H);

        gradient_t  g = gradient_t::Zero();
        hessian_t   h = hessian_t::Zero();

        double score = 0.0;
        // todo: reimplement parallelization
        for (const point_t& point_prime : points_prime)
        {
            const point_t point = t * point_prime;
            traits_t::computeGradient(map, point, J, H, score, g, h);
        }

        if (score < max_score)
        {
            lambda *= param.alpha();
            linear = linear_old;
            angular = angular_old;
            ++step_adjustments;
            continue;
        }

        if (score > max_score)
        {
            max_score = score;
            lambda = std::max(1.0, lambda / param.alpha());
            step_adjustments = 0;
        }

        /// limit H
        // cslibs_math::statistics::LimitEigenValuesByZero<DIMS>::apply(h);
        gradient_t dp = h.fullPivLu().solve(g);
        dp *= lambda;

        linear_old = linear;
        angular_old = angular;

        linear_delta = dp.template head<traits_t::LINEAR_DIMS>();
        linear += linear_delta;

        // todo: verify if we have to normalize here
        angular_delta = dp.template tail<traits_t::ANGULAR_DIMS>();
        angular += angular_delta;

        if (test_eps())
            return terminate(Termination::DELTA_EPSILON);
    }

    return terminate(Termination::MAX_ITERATIONS);
}

template<typename iterator_t, typename ndt_t, typename traits_t = MatchTraits<ndt_t>>
auto match(const ndt_t& src,
           const ndt_t& dst,
           const Parameter& param,
           const typename ndt_t::transform_t& initial_transform)
-> Result<typename ndt_t::transform_t>
{
    static constexpr int DIMS = traits_t::LINEAR_DIMS + traits_t::ANGULAR_DIMS;
    using point_t             = typename ndt_t::point_t;
    using transform_t         = typename ndt_t::transform_t;
    using result_t            = Result<transform_t>;

    using JacobianCompute = typename traits_t::Jacobian;
    using HessianCompute  = typename traits_t::Hessian;

    using linear_t      = Eigen::Matrix<double, traits_t::LINEAR_DIMS, 1>;
    using angular_t     = Eigen::Matrix<double, traits_t::ANGULAR_DIMS, 1>;
    using gradient_t    = Eigen::Matrix<double, DIMS, 1>;
    using hessian_t     = Eigen::Matrix<double, DIMS, DIMS>;

    // initialize result
    double max_score        = std::numeric_limits<double>::lowest();
    std::size_t iteration   = 0;

    linear_t  linear    = initial_transform.translation();
    angular_t angular   = initial_transform.euler();

    // initialize state
    linear_t  linear_old    = linear;
    angular_t angular_old   = angular;
    linear_t  linear_delta  = linear_t::Constant(std::numeric_limits<double>::max());
    angular_t angular_delta = angular_t::Constant(std::numeric_limits<double>::max());

    double lambda = 1.0;
    std::size_t step_adjustments = 0;

    // termination criteria
    const auto test_eps = [&]()
    {
        return (linear_delta.array().abs() < param.translationEpsilon()).all()
                && (angular_delta.array().abs() < param.rotationEpsilon()).all();
    };

    const auto test_readjustments = [&]()
    {
        return step_adjustments > 0 && step_adjustments > param.maxStepReadjustments();
    };

    // termination
    const auto terminate = [&](Termination reason)
    {
        return result_t{
            max_score,
                    iteration,
                    traits_t::makeTransform(linear, angular) * initial_transform,
                    reason };
    };


    // iterations
    for (iteration = 0; iteration < param.maxIterations(); ++iteration)
    {
        if (test_readjustments())
            return terminate(Termination::MAX_STEP_READJUSTMENTS);

        const auto t = traits_t::makeTransform(linear, angular);

        JacobianCompute J;
        JacobianCompute::get(angular, J);
        HessianCompute H;
        HessianCompute::get(angular, H);

        gradient_t  g = gradient_t::Zero();
        hessian_t   h = hessian_t::Zero();

        double score = 0.0;
        auto process_bundle = [&dst, &J, &H, &t, &score, &g, &h](const typename ndt_t::index_t &, const typename ndt_t::distribution_bundle_t &b)
        {
            traits_t::computeGradient(dst, b, J, H, t, score, g, h);
        };
        src.traverse(process_bundle);

        if (score < max_score)
        {
            lambda *= param.alpha();
            linear = linear_old;
            angular = angular_old;
            ++step_adjustments;
            continue;
        }

        if (score > max_score)
        {
            max_score = score;
            lambda = std::max(1.0, lambda / param.alpha());
            step_adjustments = 0;
        }

        /// limit H
        // cslibs_math::statistics::LimitEigenValuesByZero<DIMS>::apply(h);
        gradient_t dp = h.fullPivLu().solve(g);
        dp *= lambda;

        linear_old = linear;
        angular_old = angular;

        linear_delta = dp.template head<traits_t::LINEAR_DIMS>();
        linear += linear_delta;

        // todo: verify if we have to normalize here
        angular_delta = dp.template tail<traits_t::ANGULAR_DIMS>();
        angular += angular_delta;

        if (test_eps())
            return terminate(Termination::DELTA_EPSILON);
    }

    return terminate(Termination::MAX_ITERATIONS);
}


}
}
