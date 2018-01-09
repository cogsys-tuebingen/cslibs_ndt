#ifndef CSLIBS_NDT_INDEXED_HPP
#define CSLIBS_NDT_INDEXED_HPP

namespace cslibs_ndt {
template <template <std::size_t> class T, std::size_t Size, std::size_t Dim>
struct Indexed {
    std::array<int, Dim> index_;
    T<Size>              data_;

    Indexed() = default;

    Indexed(const std::array<int, Dim> & index,
            const T<Size>              & data) :
        index_(index),
        data_(data)
    {
    }
};
}

#endif // CSLIBS_NDT_INDEXED_HPP
