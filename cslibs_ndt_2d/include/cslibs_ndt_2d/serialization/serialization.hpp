#ifndef CSLIBS_NDT_2D_SERIALIZATION_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_HPP

#include <cslibs_ndt/map/map.hpp>
#include <cslibs_ndt/serialization/map.hpp>

namespace cslibs_ndt_2d {
namespace serialization {

template <cslibs_ndt::map::tags::option option_t,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
inline bool saveBinary(const cslibs_ndt::map::Map<option_t,2,data_t,T,backend_t,dynamic_backend_t> &map,
                       const std::string &path)
{
    return cslibs_ndt::serialization::saveBinary<option_t,2,data_t,T,backend_t,dynamic_backend_t>(map,path);
}

template <cslibs_ndt::map::tags::option option_t,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt::map::Map<option_t,2,data_t,T,backend_t,dynamic_backend_t>::Ptr& map)
{
    return cslibs_ndt::serialization::loadBinary<option_t,2,data_t,T,backend_t,dynamic_backend_t>(path,map);
}

// TODO: load with only one template arg?

}
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_HPP
