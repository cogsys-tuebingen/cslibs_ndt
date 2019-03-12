#ifndef CSLIBS_NDT_2D_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_DYNAMIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/map/impl/gridmap.hpp>

namespace cslibs_ndt_2d {
namespace dynamic_maps {

template <typename T>
using Gridmap = cslibs_ndt::map::Gridmap<cslibs_ndt::map::tags::dynamic_map,2,T>;

}
}

#endif // CSLIBS_NDT_2D_DYNAMIC_MAPS_GRIDMAP_HPP
