#ifndef CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
#define CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP

namespace cslibs_ndt_2d {
namespace matching {
namespace ceres {

enum class Flag { DIRECT, INTERPOLATION };

template <typename ndt_t, Flag flag>
class ScanMatchCostFunctor;

}
}
}

// partial specializations
#include <cslibs_ndt_2d/matching/ceres/impl/gridmap_cost_functor.hpp>
#include <cslibs_ndt_2d/matching/ceres/impl/occupancy_gridmap_cost_functor.hpp>

#endif // CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
