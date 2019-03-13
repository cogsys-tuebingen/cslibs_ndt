#ifndef CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt_2d/static_maps/gridmap.hpp>
#include <cslibs_ndt/serialization/map.hpp>

namespace cslibs_ndt_2d {
namespace static_maps {

template <typename T>
inline bool saveBinary(const typename cslibs_ndt_2d::static_maps::Gridmap<T>::Ptr &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<cslibs_ndt::map::tags::static_map,2,cslibs_ndt::Distribution,T>(map, path);
}

template <typename T>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt_2d::static_maps::Gridmap<T>::Ptr &map)
{
    return cslibs_ndt::serialization::loadBinary<cslibs_ndt::map::tags::static_map,2,cslibs_ndt::Distribution,T>(path, map);
}

}
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
