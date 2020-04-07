#ifndef CSLIBS_NDT_UTILITY_TO_POINT_HPP
#define CSLIBS_NDT_UTILITY_TO_POINT_HPP

#include <cslibs_math/linear/vector.hpp>

namespace cslibs_ndt {
namespace utility {

template <typename T, std::size_t DD>
using vector_t = cslibs_math::linear::Vector<T,DD>;

template <typename point_t, typename T, std::size_t DD, typename std::size_t... counter>
static inline point_t to_point(const vector_t<T,DD> &p, utility::integer_sequence<std::size_t,counter...>)
{
    auto at = [&p](const std::size_t &c) {
        return (c >= DD) ? T(0) : p(c);
    };
    return point_t(at(counter)...);
}

template <typename point_t, typename T, std::size_t DD>
static inline point_t to_point(const vector_t<T,DD> &p)
{
    return to_point<point_t,T,DD>(p, utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
}

template <typename point_t, typename Fn, typename std::size_t... counter>
static inline point_t to_point(const Fn &f, utility::integer_sequence<std::size_t,counter...>)
{
    auto at = [&f](const std::size_t &c) {
        return f(c);
    };
    return point_t(at(counter)...);
}

template <typename point_t, typename Fn>
static inline point_t to_point(const Fn &f)
{
    return to_point<point_t,Fn>(f, utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
}

template <std::size_t D, typename Fn, typename std::size_t... counter>
static inline std::array<int,D> to_index(const Fn &f, utility::integer_sequence<std::size_t,counter...>)
{
    auto at = [&f](const std::size_t &c) {
        return f(c);
    };
    return std::array<int,D>{{at(counter)...}};
}

template <std::size_t D, typename Fn>
static inline std::array<int,D> to_index(const Fn &f)
{
    return to_index<D,Fn>(f, utility::make_integer_sequence<std::size_t, D>{});
}

}
}

#endif // CSLIBS_NDT_UTILITY_TO_POINT_HPP
