#ifndef CSLIBS_NDT_MAP_TRAITS_HPP
#define CSLIBS_NDT_MAP_TRAITS_HPP

#include <cslibs_ndt/map/integer_sequence.hpp>

#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_math_2d/linear/pose.hpp>
#include <cslibs_math_2d/linear/pointcloud.hpp>
#include <cslibs_math_2d/algorithms/simple_iterator.hpp>

#include <cslibs_math_3d/linear/pose.hpp>
#include <cslibs_math_3d/linear/pointcloud.hpp>
#include <cslibs_math_3d/algorithms/simple_iterator.hpp>

#include <cslibs_indexed_storage/backends.hpp>
namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
namespace map {

namespace tags {
struct option {};

struct static_map : option {
    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using default_backend_t         = cis::backend::array::Array<data_interface_t_, index_interface_t_, options_ts_...>;
    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using default_dynamic_backend_t = cis::backend::kdtree::KDTree<data_interface_t_, index_interface_t_, options_ts_...>;
};
struct dynamic_map : option {
    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using default_backend_t         = cis::backend::kdtree::KDTree<data_interface_t_, index_interface_t_, options_ts_...>;
    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using default_dynamic_backend_t = cis::backend::kdtree::KDTree<data_interface_t_, index_interface_t_, options_ts_...>;
};
}

template <std::size_t Dim, typename T>
struct traits {};

template <typename T>
struct traits<2,T>
{
    using pose_t                = cslibs_math_2d::Pose2d<T>;
    using transform_t           = cslibs_math_2d::Transform2d<T>;
    using point_t               = cslibs_math_2d::Point2d<T>;
    using pointcloud_t          = cslibs_math_2d::Pointcloud2d<T>;
    using default_iterator_t    = cslibs_math_2d::algorithms::SimpleIterator<T>;
};

template <typename T>
struct traits<3,T>
{
    using pose_t                = cslibs_math_3d::Pose3d<T>;
    using transform_t           = cslibs_math_3d::Transform3d<T>;
    using point_t               = cslibs_math_3d::Point3d<T>;
    using pointcloud_t          = cslibs_math_3d::Pointcloud3d<T>;
    using default_iterator_t    = cslibs_math_3d::algorithms::SimpleIterator<T>;
};

namespace detail {
inline constexpr std::size_t two_pow(const std::size_t Dim)
{
    return (Dim == 0) ? 1 : (2 * two_pow(Dim-1));
}

template <typename index_t>
inline constexpr int get_index(const index_t bi, const std::size_t counter, const std::size_t dimension)
{
    return ((counter / two_pow(dimension)) % 2 == 0) ?
              cslibs_math::common::div(bi[dimension],2) :
              (cslibs_math::common::div(bi[dimension],2) + cslibs_math::common::mod(bi[dimension],2));
}
template <typename index_t, std::size_t... dimension>
inline constexpr index_t generate_index(const index_t bi, const std::size_t counter, integer_sequence<std::size_t,dimension...>)
{
    return {get_index<index_t>(bi,counter,dimension)...};
}
template <typename index_t,
          typename integers = make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>>
inline constexpr index_t generate_index(const index_t bi, const std::size_t counter)
{
    return generate_index<index_t>(bi,counter,integers{});
}
template <typename list_t, typename index_t, std::size_t... counter>
inline constexpr list_t generate_indices_helper(const index_t bi, integer_sequence<std::size_t,counter...>)
{
    return {generate_index<index_t>(bi,counter)...};
}
template <typename list_t,
          std::size_t Dim,
          typename index_t = std::array<int,Dim>,
          typename integers = make_integer_sequence<std::size_t,std::tuple_size<list_t>::value>>
inline constexpr list_t generate_indices(const index_t bi)
{
    return generate_indices_helper<list_t,index_t>(bi,integers{});
}


inline bool merge_and()
{
    return true;
}
template <typename... V>
inline bool merge_and(const bool& v, const V&... values)
{
    return v && merge_and(values...);
}
}

}
}

#endif // CSLIBS_NDT_MAP_TRAITS_HPP
