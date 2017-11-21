#ifndef MASK_HPP
#define MASK_HPP

#include <vector>
#include <array>
#include <iostream>
#include <assert.h>

namespace ndt {
namespace grid {
template<std::size_t Iteration>
struct pow2 {
    static const std::size_t value = 2 * pow2<Iteration - 1>::value;
};
template<>
struct pow2<0> {
    static const std::size_t value = 1;
};

template<std::size_t Iteration, typename Mask>
struct fill {
    static constexpr void assign(Mask &mask)
    {
        const std::size_t off = mask.cols - Iteration;
        const std::size_t div = pow2<Iteration - 1>::value;
        for(std::size_t i = 0 ; i < mask.rows; ++i) {
            mask[i * mask.cols + off] = (i / div) % 2;
//            if(((i / div) % 2) == 0)
//                mask[i * mask.cols + off] = -1;
//            else
//                mask[i * mask.cols + off] = 1;
        }
        fill<Iteration - 1, Mask>::assign(mask);
    }

};

template<typename Mask>
struct fill<0, Mask>
{
    static constexpr void assign(Mask &mask)
    {
        return;
    }
};

/**
 * @brief The Mask struct can be used for binary couting depending on the dimension given to it.
 */
template<std::size_t Dim>
struct Mask {
    static const std::size_t rows = pow2<Dim>::value;
    static const std::size_t cols = Dim;
    std::array<int, rows * cols> mask;
    std::size_t                  pos[rows];

    inline int & operator [] (const std::size_t i)
    {
        return mask[i];
    }

    inline const int & operator [] (const std::size_t i) const
    {
        return mask[i];
    }

    Mask()
    {
        /// fill the mask
        fill<Dim, Mask<Dim>>::assign(*this);
        for(std::size_t i = 0 ; i < rows ; ++i) {
            pos[i] = i * cols;
        }
    }
};
}
}
#endif // MASK_HPP
