#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP

namespace cslibs_ndt {
namespace matching {
namespace ceres {

enum class Flag { DIRECT, INTERPOLATION };

template <typename ndt_t, Flag flag>
class ScanMatchCostFunctor;

}
}
}

// partial specializations
#include <cslibs_ndt/matching/ceres/map/impl/gridmap_cost_functor.hpp>
#include <cslibs_ndt/matching/ceres/map/impl/occupancy_gridmap_cost_functor.hpp>

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
