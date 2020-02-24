#ifndef CSLIBS_NDT_UTILITY_BILINEAR_INTERPOLATION_HPP
#define CSLIBS_NDT_UTILITY_BILINEAR_INTERPOLATION_HPP

#include <cslibs_ndt/utility/integer_sequence.hpp>
#include <cslibs_ndt/utility/bilinear_interpolation.hpp>

#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <type_traits>
#include <numeric>

namespace cslibs_ndt {
namespace utility {

template <typename T, typename std::size_t D, typename point_t, typename std::size_t... counter>
static inline std::array<T,D> get_bilinear_interpolation_weights(const std::array<int,D>& i,
                                                                 const point_t& p,
                                                                 const T &res_inv,
                                                                 utility::integer_sequence<std::size_t,counter...>)
{
    auto at = [&i,&p,&res_inv](const std::size_t &c) {
        return (c >= D) ? T(0) :
                          (i[c] & 1) ? (((p(c) * res_inv) - static_cast<T>(i[c]))) :
                                       (T(1.) - ((p(c) * res_inv) - static_cast<T>(i[c])));
    };
    return std::array<T,D>{at(counter)...};
}

template <typename T, typename std::size_t D, typename point_t>
static inline std::array<T,D> get_bilinear_interpolation_weights(const std::array<int,D>& i,
                                                                 const point_t& p,
                                                                 const T& res_inv)
{
    return get_bilinear_interpolation_weights<T,D,point_t>(
                i,p,res_inv,utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
}

template <typename T>
static inline T to_bilinear_interpolation_weight_inner(const T& w, const std::size_t& counter, const std::size_t& dimension)
{
    return ((counter / two_pow(dimension)) & 1) ? w : (T(1.) - w);
}


template <typename T, std::size_t D>
static inline T to_bilinear_interpolation_weight(const std::array<T,D>& w, const std::size_t& counter, const std::size_t& dimension)
{
    return to_bilinear_interpolation_weight_inner(w[dimension], counter, dimension);
}

template <typename T, std::size_t D, std::size_t... dimension>
static inline std::array<T,D> to_bilinear_interpolation_weights(const std::array<T,D>& w, const std::size_t& counter, integer_sequence<std::size_t,dimension...>)
{
    return {to_bilinear_interpolation_weight<T,D>(w,counter,dimension)...};
}

template <typename T, std::size_t D,
          typename integers = make_integer_sequence<std::size_t,D>>
static inline std::array<T,D> to_bilinear_interpolation_weights(const std::array<T,D>& w, const std::size_t& counter)
{
    return to_bilinear_interpolation_weights<T,D>(w,counter,integers{});
}

template <typename T, std::size_t D>
static inline T to_bilinear_interpolation_weight(const std::array<T,D>& w, const std::size_t& counter)
{
    const std::array<T,D>& weights = to_bilinear_interpolation_weights(w,counter);
    return std::accumulate(weights.begin(), weights.end(), T(1), std::multiplies<T>());
}

}
}

#endif // CSLIBS_NDT_UTILITY_BILINEAR_INTERPOLATION_HPP
