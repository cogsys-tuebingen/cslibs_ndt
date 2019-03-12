#ifndef CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/map/impl/gridmap.hpp>

namespace cslibs_ndt_3d {
namespace static_maps {

template <typename T>
using Gridmap = cslibs_ndt::map::Gridmap<cslibs_ndt::map::tags::static_map,3,T>;

}
}

#endif // CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP
