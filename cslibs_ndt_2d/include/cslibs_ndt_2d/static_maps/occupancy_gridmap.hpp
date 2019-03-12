#ifndef CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/map/impl/occupancy_gridmap.hpp>

namespace cslibs_ndt_2d {
namespace static_maps {

template <typename T>
using OccupancyGridmap = cslibs_ndt::map::OccupancyGridmap<cslibs_ndt::map::tags::static_map,2,T>;

}
}

#endif // CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
