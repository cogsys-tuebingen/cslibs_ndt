#pragma once

#include <cslibs_ndt/matching/match_traits.hpp>
#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/matching/jacobian.hpp>
#include <cslibs_ndt_3d/matching/hessian.hpp>

namespace cslibs_ndt {
namespace matching {

template<>
struct MatchTraits<cslibs_ndt_3d::dynamic_maps::Gridmap>
{
    static constexpr int LINEAR_DIMS  = 3;
    static constexpr int ANGULAR_DIMS = 3;
    using Jacobian = cslibs_ndt_3d::matching::Jacobian;
    using Hessian  = cslibs_ndt_3d::matching::Hessian;
};

}
}
