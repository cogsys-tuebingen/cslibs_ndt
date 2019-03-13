#ifndef CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt/serialization/generic_map.hpp>

namespace cslibs_ndt_2d {
namespace dynamic_maps {

template <typename T>
inline bool saveBinary(const typename cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>::Ptr &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<cslibs_ndt::map::tags::dynamic_map,2,cslibs_ndt::OccupancyDistribution,T>(map, path);
}

template <typename T>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>::Ptr &map)
{
    using target_ptr_t = typename cslibs_ndt::map::GenericMap<cslibs_ndt::map::tags::dynamic_map,2,cslibs_ndt::OccupancyDistribution,T,
    cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::template default_backend_t,
    cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::template default_dynamic_backend_t>::Ptr;
    target_ptr_t &m = reinterpret_cast<target_ptr_t&>(map);
    return cslibs_ndt::serialization::loadBinary<cslibs_ndt::map::tags::dynamic_map,2,cslibs_ndt::OccupancyDistribution,T>(path, m);
}

}
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
