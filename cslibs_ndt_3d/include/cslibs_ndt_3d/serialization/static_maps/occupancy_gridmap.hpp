#ifndef CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt_3d/static_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt/serialization/generic_map.hpp>

namespace cslibs_ndt_3d {
namespace static_maps {

template <typename T>
inline bool saveBinary(const typename cslibs_ndt_3d::static_maps::OccupancyGridmap<T>::Ptr &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<cslibs_ndt::map::tags::static_map,3,cslibs_ndt::OccupancyDistribution,T>(map, path);
}

template <typename T>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt_3d::static_maps::OccupancyGridmap<T>::Ptr &map)
{
    using target_ptr_t = typename cslibs_ndt::map::GenericMap<cslibs_ndt::map::tags::static_map,3,cslibs_ndt::OccupancyDistribution,T,
    cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::static_map>::template default_backend_t,
    cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::static_map>::template default_dynamic_backend_t>::Ptr;
    target_ptr_t &m = reinterpret_cast<target_ptr_t&>(map);
    return cslibs_ndt::serialization::loadBinary<cslibs_ndt::map::tags::static_map,3,cslibs_ndt::OccupancyDistribution,T>(path, m);
}

}
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
