#ifndef CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt/serialization/map.hpp>

namespace cslibs_ndt_3d {
namespace dynamic_maps {

template <typename T>
inline bool saveBinary(const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::binary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::Distribution,T>::save(*map, path);
}

template <typename T>
inline bool saveBinary(const cslibs_ndt_3d::dynamic_maps::Gridmap<T> &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::binary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::Distribution,T>::save(map, path);
}

template <typename T>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &map)
{
    return cslibs_ndt::serialization::binary<cslibs_ndt::map::tags::dynamic_map,3,cslibs_ndt::Distribution,T>::load(path, map);
}

}
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
