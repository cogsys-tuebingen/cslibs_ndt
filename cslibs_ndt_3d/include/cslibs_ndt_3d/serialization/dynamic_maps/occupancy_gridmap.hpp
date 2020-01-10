#ifndef CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt/serialization/map.hpp>

namespace cslibs_ndt_3d {
namespace dynamic_maps {

template <typename T>
inline bool saveBinary(const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::OccupancyDistribution,T>(*map, path);
}

template <typename T>
inline bool saveBinary(const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T> &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::OccupancyDistribution,T>(map, path);
}

template <typename T>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &map)
{
    return cslibs_ndt::serialization::loadBinary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::OccupancyDistribution,T>(path, map);
}

}
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
