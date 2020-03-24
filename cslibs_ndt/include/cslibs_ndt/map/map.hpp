#ifndef CSLIBS_NDT_MAP_MAP_HPP
#define CSLIBS_NDT_MAP_MAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>

namespace cslibs_ndt {
namespace map {

template <tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = tags::default_types<option_t>::template default_backend_t>
class EIGEN_ALIGN16 Map : public GenericMap<option_t,Dim,data_t,T,backend_t> {};

}
}

// partial specializations
#include <cslibs_ndt/map/impl/gridmap.hpp>
#include <cslibs_ndt/map/impl/occupancy_gridmap.hpp>
#include <cslibs_ndt/map/impl/weighted_occupancy_gridmap.hpp>

#endif // CSLIBS_NDT_MAP_MAP_HPP
