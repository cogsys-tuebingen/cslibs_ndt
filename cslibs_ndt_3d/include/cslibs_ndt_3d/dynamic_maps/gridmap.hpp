#ifndef CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/map/map.hpp>

namespace cslibs_ndt_3d {
namespace dynamic_maps {

template <typename T>
using Gridmap = cslibs_ndt::map::Map<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::Distribution,T>;

}
}

#endif // CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP
