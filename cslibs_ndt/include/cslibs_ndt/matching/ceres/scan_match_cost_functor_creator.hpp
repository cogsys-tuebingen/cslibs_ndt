#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_CREATOR_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_CREATOR_HPP

#include <ceres/cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/numeric_diff_cost_function.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <template <typename,typename> class child_t, typename base_t>
class ScanMatchCostFunctorCreator
{
public:
    template <typename points_t, typename ... args_t>
    static inline ::ceres::CostFunction* CreateAutoDiffCostFunction(double weight,
                                                                    points_t points,
                                                                    const args_t &...args)
    {
        using _child_t = child_t<base_t,points_t>;

        const auto count = points.size();
        return new ::ceres::AutoDiffCostFunction<_child_t, ::ceres::DYNAMIC, _child_t::N0, _child_t::N1>(
                    new _child_t(weight / std::sqrt(count), std::move(points), args...),
                    count);
    }

    template <::ceres::NumericDiffMethodType method = ::ceres::FORWARD, typename points_t, typename ... args_t>
    static inline ::ceres::CostFunction* CreateNumericDiffCostFunction(double weight,
                                                                       points_t points,
                                                                       const args_t &...args)
    {
        using _child_t = child_t<base_t,points_t>;

        const auto count = points.size();
        return new ::ceres::NumericDiffCostFunction<_child_t, method, ::ceres::DYNAMIC, _child_t::N0, _child_t::N1>(
                    new _child_t(weight / std::sqrt(count), std::move(points), args...),
                    ::ceres::DO_NOT_TAKE_OWNERSHIP,
                    count);
    }
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_CREATOR_HPP
