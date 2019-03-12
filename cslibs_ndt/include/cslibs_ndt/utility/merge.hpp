#ifndef CSLIBS_NDT_UTILITY_MERGE_HPP
#define CSLIBS_NDT_UTILITY_MERGE_HPP

namespace cslibs_ndt {
namespace utility {

template <typename Fn>
inline typename Fn::return_type merge()
{
    return Fn::default_value;
}

template <typename Fn, typename... V>
inline typename Fn::return_type merge(const typename Fn::return_type& v, const V&... values)
{
    return Fn::apply(v, merge<Fn>(values...));
}

namespace operations {
struct bool_and {
    using return_type = bool;
    static const return_type default_value = true;

    static return_type apply(const return_type& v1, const return_type& v2)
    {
        return v1 && v2;
    }
};

struct bool_or {
    using return_type = bool;
    static const return_type default_value = false;

    static return_type apply(const return_type& v1, const return_type& v2)
    {
        return v1 || v2;
    }
};
}

}
}

#endif // CSLIBS_NDT_UTILITY_MERGE_HPP
