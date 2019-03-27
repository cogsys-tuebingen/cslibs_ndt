#ifndef CSLIBS_NDT_UTILITY_FOR_EACH_HPP
#define CSLIBS_NDT_UTILITY_FOR_EACH_HPP

namespace cslibs_ndt {
namespace utility {

namespace detail {

template <std::size_t N, typename Fn>
struct for_each_helper {
    static inline void apply(Fn& function)
    {
        for_each_helper<N-1,Fn>::apply(function);
        function(N);
    }
};

template <typename Fn>
struct for_each_helper<0,Fn> {
    static inline void apply(Fn& function)
    {
        function(0);
    }
};

}

template <std::size_t N, typename Fn>
static inline void for_each(Fn &function)
{
    detail::for_each_helper<N-1,Fn>::apply(function);
}

}
}

#endif // CSLIBS_NDT_UTILITY_FOR_EACH_HPP
