#ifndef CSLIBS_NDT_UTILITY_FOR_EACH_HPP
#define CSLIBS_NDT_UTILITY_FOR_EACH_HPP

#include <type_traits>

namespace cslibs_ndt {
namespace utility {

template <std::size_t From, std::size_t To, typename Fn>
static inline typename std::enable_if<From == To, void>::type
for_each(const Fn &function)
{
}

template <std::size_t From, std::size_t To, typename Fn>
static inline typename std::enable_if<From < To, void>::type
for_each(const Fn &function)
{
    function(From);
    for_each<From+1,To,Fn>(function);
}

template <std::size_t To, typename Fn>
static inline void for_each(const Fn &function)
{
    for_each<0,To,Fn>(function);
}

}
}

#endif // CSLIBS_NDT_UTILITY_FOR_EACH_HPP
