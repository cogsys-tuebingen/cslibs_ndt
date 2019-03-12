#ifndef CSLIBS_NDT_UTILITY_CREATE_HPP
#define CSLIBS_NDT_UTILITY_CREATE_HPP

#include <memory>

namespace cslibs_ndt {
namespace utility {

template <typename storage_t, std::size_t Dim, std::size_t... counter>
inline std::array<std::shared_ptr<storage_t>,Dim> create(
        const std::array<std::shared_ptr<storage_t>,Dim>& other,
        utility::integer_sequence<std::size_t,counter...>)
{
    auto create_at = [&other](const std::size_t& c) {
        return std::shared_ptr<storage_t>(new storage_t(*other[c]));
    };
    return {create_at(counter)...};
}

template <typename storage_t, std::size_t Dim>
inline std::array<std::shared_ptr<storage_t>,Dim> create(
        const std::array<std::shared_ptr<storage_t>,Dim>& other)
{
    return create<storage_t,Dim>(other,utility::make_integer_sequence<std::size_t,Dim>{});
}


template <typename storage_t, std::size_t Dim, std::size_t... counter>
inline std::array<std::shared_ptr<storage_t>,Dim> create(
        utility::integer_sequence<std::size_t,counter...>)
{
    auto create_at = [](const std::size_t& c) {
        return std::shared_ptr<storage_t>(new storage_t);
    };
    return {create_at(counter)...};
}

template <typename storage_t, std::size_t Dim>
inline std::array<std::shared_ptr<storage_t>,Dim> create()
{
    return create<storage_t,Dim>(utility::make_integer_sequence<std::size_t,Dim>{});
}


template <typename T, std::size_t Dim, std::size_t... counter>
inline std::array<T,Dim> create(const T& value, utility::integer_sequence<std::size_t,counter...>)
{
    auto create_at = [&value](const std::size_t& c) {
        return value;
    };
    return {create_at(counter)...};
}

template <typename T, std::size_t Dim>
inline std::array<T,Dim> create(const T& value)
{
    return create<T,Dim>(value,utility::make_integer_sequence<std::size_t,Dim>{});
}

}
}

#endif // CSLIBS_NDT_UTILITY_CREATE_HPP
