#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP

namespace cslibs_ndt {
namespace matching {
namespace ceres {

enum class Flag { DIRECT, INTERPOLATION };

template <typename ndt_t, Flag flag_t>
class ScanMatchCostFunctor;

}
}
}

// partial specializations can be found in
// cslibs_ndt_2d,
// cslibs_ndt_3d

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
