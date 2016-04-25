#ifndef SGN_HPP
#define SGN_HPP

namespace ndt {
namespace math {
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}
}

#endif // SGN_HPP
