#ifndef CSLIBS_NDT_MATCHING_CERES_LOCAL_PARAMETERIZATION_HPP
#define CSLIBS_NDT_MATCHING_CERES_LOCAL_PARAMETERIZATION_HPP

#include <ceres/autodiff_local_parameterization.h>
#include <ceres/rotation.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

struct YawOnlyQuaternionPlus
{
    static inline ::ceres::LocalParameterization* CreateAutoDiff()
    {
        return new ::ceres::AutoDiffLocalParameterization<YawOnlyQuaternionPlus, 4, 1>();
    }

    template<typename T>
    inline bool operator()(const T* x, const T* delta, T* x_plus_delta) const
    {
        auto clamp = [](const T& value, const T& min, const T& max) {
            if (value > max) return max;
            if (value < min) return min;
            return value;
        };

        const T clamped_delta = clamp(delta[0], T(-0.5), T(0.5));
        T q_delta[4];
        q_delta[0] = ::ceres::sqrt(1. - clamped_delta * clamped_delta);
        q_delta[1] = T(0.);
        q_delta[2] = T(0.);
        q_delta[3] = clamped_delta;
        ::ceres::QuaternionProduct(q_delta, x, x_plus_delta);
        return true;
    }
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_LOCAL_PARAMETERIZATION_HPP
