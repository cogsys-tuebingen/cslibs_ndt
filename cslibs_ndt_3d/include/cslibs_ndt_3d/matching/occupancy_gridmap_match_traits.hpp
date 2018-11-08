#pragma once

#include <cslibs_ndt/matching/match_traits.hpp>
#include <cslibs_ndt/matching/occupancy_parameter.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_3d/matching/jacobian.hpp>
#include <cslibs_ndt_3d/matching/hessian.hpp>

namespace cslibs_ndt {
namespace matching {

template<typename MapT> struct IsOccupancyGridmap : std::false_type {};
template<> struct IsOccupancyGridmap<cslibs_ndt_3d::dynamic_maps::OccupancyGridmap> : std::true_type {};
template<> struct IsOccupancyGridmap<cslibs_ndt_3d::static_maps::OccupancyGridmap> : std::true_type {};

template<typename MapT>
struct MatchTraits<MapT, typename std::enable_if<IsOccupancyGridmap<MapT>::value>::type>
{
    static constexpr int LINEAR_DIMS  = 3;
    static constexpr int ANGULAR_DIMS = 3;
    using Jacobian  = cslibs_ndt_3d::matching::Jacobian;
    using Hessian   = cslibs_ndt_3d::matching::Hessian;

    using gradient_t = Eigen::Matrix<double, 6, 1>;
    using hessian_t  = Eigen::Matrix<double, 6, 6>;

    using point_t = cslibs_math_3d::Point3d;
    using transform_t = cslibs_math_3d::Transform3d;
    using parameter_t = cslibs_ndt::matching::OccupancyParameter;

    static transform_t makeTransform(const Eigen::Vector3d& linear,
                                     const Eigen::Vector3d& angular)
    {
        return transform_t{
                linear.x(), linear.y(), linear.z(),
                angular.x(), angular.y(), angular.z()};
    }

    // todo: deduplicate code, make model configureable...
    static void computeGradient(const MapT& map,
                                const point_t& point,
                                const Jacobian& J,
                                const Hessian& H,
                                const parameter_t& param,
                                double& score,
                                gradient_t& g,
                                hessian_t& h)
    {
        static constexpr double d1 = 0.95;
        static constexpr double d2 = 1 - d1;

        auto* bundle = map.getDistributionBundle(point);
        if (!bundle)
            return;

        // check occupancy value
        if (param.occupancyThreshold() > 0.0)
        {
            double occupancy = 0.0;
            for (auto* distribution_wrapper : *bundle)
                occupancy += distribution_wrapper->getOccupancy(param.inverseModel());
            occupancy /= 8.0;

            if (occupancy < param.occupancyThreshold())
                return;
        }

        for (auto* distribution_wrapper : *bundle)
        {
            auto& d = distribution_wrapper->getDistribution();
            if (!d || d->getN() < 4)
                continue;

            const auto info   = d->getInformationMatrix();
            const auto q      = (point.data() - d->getMean()).eval();
            const auto q_info = (q.transpose() * info).eval();
            const auto p_occ  = distribution_wrapper->getOccupancy(param.inverseModel()); // no recompute: this uses a cached value
            const auto e      = -0.5 * double(q_info * q) * (d2 * (1 - p_occ));
            const auto s      = d1 * p_occ * std::exp(e);
            if (!std::isnormal(s) || s <= 1e-5)
                continue;

            // this part should be vectorized, may also remove common factors...
            for (std::size_t i = 0; i < LINEAR_DIMS + ANGULAR_DIMS; ++i)
            {
                const auto J_iq = J.get(i, q);
                const auto J_info = (info * J_iq).eval();

                g(i) += s * q_info * J_iq;

                for (std::size_t j = 0; j < LINEAR_DIMS + ANGULAR_DIMS; ++j)
                {
                    h(i, j) -= s * q_info * H.get(i, j, q) +
                               s * static_cast<double>((J.get(j, q).transpose()).eval() * J_info) -
                               s * (q_info * J_iq).eval() * (-q_info * J.get(j, q)).eval();
                }
            }

            score += s;
        }
    }
};

}
}
