#ifndef CSLIBS_NDT_UTILITY_BINARY_INDICES_HPP
#define CSLIBS_NDT_UTILITY_BINARY_INDICES_HPP

#include <cslibs_ndt/utility/integer_sequence.hpp>

#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

namespace cslibs_ndt {
namespace utility {

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

}
}

#endif // CSLIBS_NDT_UTILITY_BINARY_INDICES_HPP
