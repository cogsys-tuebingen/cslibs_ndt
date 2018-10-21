#pragma once

#include <cslibs_ndt_3d/matching/gridmap_match_traits.hpp>
#include <cslibs_ndt/matching/match_traits.hpp>

namespace cslibs_ndt {
namespace matching {

struct Parameter
{
    cslibs_math_3d::Transform3d transform;
    double trans_eps;
    double rot_eps;
    std::size_t max_step_readjustments;
    std::size_t max_iterations;
    double alpha;
};

struct Result
{
    enum class Termination { ITER, EPS, READJUST };

    double score;
    std::size_t iterations;
    cslibs_math_3d::Transform3d transform;
    Termination termination;
};

template<typename iterator_t>//, typename ndt_t, typename traits_t = MatchTraits<ndt_t>>
Result match(const iterator_t& points_begin,
             const iterator_t& points_end,
             const cslibs_ndt_3d::dynamic_maps::Gridmap& map,
             const Parameter& param)
{
    using ndt_t     = cslibs_ndt_3d::dynamic_maps::Gridmap;
    using traits_t  = MatchTraits<cslibs_ndt_3d::dynamic_maps::Gridmap>;

    static constexpr int DIMS = traits_t::LINEAR_DIMS + traits_t::ANGULAR_DIMS;
    using point_t             = typename ndt_t::point_t;
    using point_transposed_t  = Eigen::Matrix<double, 1, 3>;
    using transform_t         = typename ndt_t::transform_t;
    using bundle_t            = typename ndt_t::distribution_bundle_t;
    using ndt_distribution_t  = typename ndt_t::distribution_t;
    using distribution_t      = typename ndt_distribution_t::distribution_t;
    using distributions_t     = std::array<distribution_t, 8>;

    using JacobianCompute = typename traits_t::Jacobian;
    using HessianCompute  = typename traits_t::Hessian;

    using linear_t      = Eigen::Matrix<double, traits_t::LINEAR_DIMS, 1>;
    using angular_t     = Eigen::Matrix<double, traits_t::ANGULAR_DIMS, 1>;
    using gradient_t    = Eigen::Matrix<double, DIMS, 1>;
    using hessian_t     = Eigen::Matrix<double, DIMS, DIMS>;

    // pre transform points, should be externalized or made completely optional...
    std::vector<point_t> points_prime;
    points_prime.reserve(std::distance(points_begin, points_end));
    std::transform(points_begin, points_end, std::back_inserter(points_prime),
            [&](const point_t& point) { return param.transform * point; });

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
        return (linear_delta.array().abs() < param.trans_eps).all()
               && (angular_delta.array().abs() < param.rot_eps).all();
    };

    const auto test_readjustments = [&]()
    {
        return step_adjustments > 0 && step_adjustments > param.max_step_readjustments;
    };

    // termination
    const auto terminate = [&](Result::Termination reason)
    {
        return Result{
            max_score,
            iteration,
            cslibs_math_3d::Transform3d{
                linear.x(), linear.y(), linear.z(),
                angular.x(), angular.y(), angular.z()} * param.transform,
            reason };
    };

    // iterations
    for (iteration = 0; iteration < param.max_iterations; ++iteration)
    {
        if (test_readjustments())
            return terminate(Result::Termination::READJUST);

        const transform_t t(linear.x(), linear.y(), linear.z(), angular.x(), angular.y(), angular.z());

        JacobianCompute J;
        JacobianCompute::get(angular, J);
        HessianCompute H;
        HessianCompute::get(angular, H);

        gradient_t  g = gradient_t::Zero();
        hessian_t   h = hessian_t::Zero();

        double score = 0.0;
        for (const point_t& point_prime : points_prime)
        {
            const point_t point = t * point_prime;
            auto* bundle = map.getDistributionBundle(point);
            if (!bundle)
                continue;

            for (auto* distribution_wrapper : *bundle)
            {
                auto& d = distribution_wrapper->data();
                if (d.getN() < 3)
                    continue;

                const auto info   = d.getInformationMatrix();
                const auto q      = (point.data() - d.getMean()).eval();
                const auto q_info = (q.transpose() * info).eval();
                const auto e      = -0.5 * double(q_info * q);
                const auto s      = std::exp(e);

                for (std::size_t i = 0; i < DIMS; ++i)
                {
                    const auto J_iq = J.get(i, q);
                    const auto J_info = (J_iq.transpose() * info).eval();

                    g(i) += s * q_info * J_iq;

                    for (std::size_t j = 0; j < DIMS; ++j)
                    {
                        h(i, j) += s * q_info * H.get(i, j, q) +
                                s * static_cast<double>(J_info * J.get(j, q));
                    }
                }

                score += s;
            }
        }

        if (score < max_score)
        {
            lambda *= param.alpha;
            linear = linear_old;
            angular = angular_old;
            ++step_adjustments;
            continue;
        }

        if (score > max_score)
        {
            max_score = score;
            lambda /= param.alpha;
            step_adjustments = 0;
        }

        cslibs_math::statistics::LimitEigenValues<DIMS, 3>::apply(h);
        gradient_t dp = -h.fullPivLu().solve(g);
        dp *= lambda;

        linear_old = linear;
        angular_old = angular;

        linear_delta = dp.template head<traits_t::LINEAR_DIMS>();
        linear += linear_delta;
        
        angular_delta = dp.template tail<traits_t::ANGULAR_DIMS>();
        angular += angular_delta;

        if (test_eps())
            return terminate(Result::Termination::EPS);
    }

    return terminate(Result::Termination::ITER);
}

}
}
