#ifndef ANGLE_HPP
#define ANGLE_HPP
#include <cmath>

namespace ndt {
namespace math {
#define _2_M_PI 2.0 * M_PI

inline double wrapAngle(const double angle )
{
    return angle - _2_M_PI * floor( angle / _2_M_PI );
}
}
}

#endif // ANGLE_HPP
