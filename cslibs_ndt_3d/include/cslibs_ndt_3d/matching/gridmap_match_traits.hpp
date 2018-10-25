#pragma once

#include <cslibs_ndt/matching/match_traits.hpp>
#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/gridmap.hpp>
#include <cslibs_ndt_3d/matching/jacobian.hpp>
#include <cslibs_ndt_3d/matching/hessian.hpp>

namespace cslibs_ndt {
namespace matching {

template<typename MapT> struct IsGridmap : std::false_type {};
template<> struct IsGridmap<cslibs_ndt_3d::dynamic_maps::Gridmap> : std::true_type {};
template<> struct IsGridmap<cslibs_ndt_3d::static_maps::Gridmap> : std::true_type {};

template<typename MapT>
struct MatchTraits<MapT, typename std::enable_if<IsGridmap<MapT>::value>::type>
{
    static constexpr int LINEAR_DIMS  = 3;
    static constexpr int ANGULAR_DIMS = 3;
    using Jacobian  = cslibs_ndt_3d::matching::Jacobian;
    using Hessian   = cslibs_ndt_3d::matching::Hessian;

    using gradient_t = Eigen::Matrix<double, 6, 1>;
    using hessian_t  = Eigen::Matrix<double, 6, 6>;

    using point_t = cslibs_math_3d::Point3d;
    using transform_t = cslibs_math_3d::Transform3d;

    static transform_t makeTransform(const Eigen::Vector3d& linear,
                                     const Eigen::Vector3d& angular)
    {
        return transform_t{
                linear.x(), linear.y(), linear.z(),
                angular.x(), angular.y(), angular.z()};
    }

    static void computeGradient(const MapT& map,
                                const point_t& point,
                                const Jacobian& J,
                                const Hessian& H,
                                double& score,
                                gradient_t& g,
                                hessian_t& h)
    {
        auto* bundle = map.getDistributionBundle(point);
        if (!bundle)
            return;

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

            // this part should be vectorized, may also remove common factors...
            for (std::size_t i = 0; i < LINEAR_DIMS + ANGULAR_DIMS; ++i)
            {
                const auto J_iq = J.get(i, q);
                const auto J_info = (J_iq.transpose() * info).eval();

                g(i) += s * q_info * J_iq;

                for (std::size_t j = 0; j < LINEAR_DIMS + ANGULAR_DIMS; ++j)
                {
                    h(i, j) += s * q_info * H.get(i, j, q) +
                               s * static_cast<double>(J_info * J.get(j, q));
                }
            }

            score += s;
        }
    }
};

}
}
