#ifndef CSLIBS_NDT_MAP_INTEGER_SEQUENCE_HPP
#define CSLIBS_NDT_MAP_INTEGER_SEQUENCE_HPP

namespace cslibs_ndt {
namespace map {
namespace detail {

template <typename T, T... I>
struct integer_sequence
{
    using type = integer_sequence;
    using value_type = T;
    static constexpr std::size_t size() noexcept { return sizeof...(I); }
};

template<class Sequence1, class Sequence2>
struct merge_and_renumber;

template<typename T, T... I1, T... I2>
struct merge_and_renumber<integer_sequence<T, I1...>, integer_sequence<T, I2...>>
        : integer_sequence<T, I1..., (sizeof...(I1) + I2)...>
{};

template <typename T, std::size_t N>
struct integer_sequence_maker
        : detail::merge_and_renumber<typename integer_sequence_maker<T, N/2>::type,
                typename integer_sequence_maker<T, N - N/2>::type>
{};

template<typename T> struct integer_sequence_maker<T,0> : integer_sequence<T> { };
template<typename T> struct integer_sequence_maker<T,1> : integer_sequence<T,0> { };

template<typename T, std::size_t N>
using make_integer_sequence = typename integer_sequence_maker<T,N>::type;
}
}
}

#endif // CSLIBS_NDT_MAP_INTEGER_SEQUENCE_HPP
