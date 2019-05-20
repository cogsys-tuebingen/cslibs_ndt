#ifndef CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_QUATERNION_3D_HPP
#define CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_QUATERNION_3D_HPP

#include <cslibs_math_3d/linear/quaternion.hpp>

#include <ceres/cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/jet.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

class RotationCostFunctor3dRPY
{
public:
    static ::ceres::CostFunction* CreateAutoDiffCostFunction(double weight,
                                                             const double& roll,
                                                             const double& pitch,
                                                             const double& yaw)
    {
        return new ::ceres::AutoDiffCostFunction<RotationCostFunctor3dRPY, 3, 3>(
                new RotationCostFunctor3dRPY(weight, roll, pitch, yaw)
        );
    }

    template<typename T>
    bool operator()(const T* const rotation_rpy, T* residual) const
    {
        for (std::size_t i=0; i<3; ++i) {
            T delta;
            T r0(rpy_[i]);
            AngleDifference(&r0, rotation_rpy[i], &delta);

            residual[i] = weight_ * delta;
        }
        return true;
    }

private:
    template<typename T>
    static void AngleDifference(const T* const r0, const T* const r, T* const d)
    {
        static const double _2_M_PI = 2.0 * M_PI;
        auto norm = [](const T* const r) {
            return ::ceres::atan2(::ceres::sin(r[0]), ::ceres::cos(r[0]));
        };

        const auto a = norm(r0);
        const auto b = norm(r);

        const auto d1 = a - b;
        const auto d2 = (_2_M_PI - ::ceres::abs(d1)) * (d1 > T(0) ? T(-1) : T(1));
        *d = ::ceres::abs(d1) < ::ceres::abs(d2) ? d1 : d2;
    }

    explicit RotationCostFunctor3dRPY(double weight,
                                      const double& roll,
                                      const double& pitch,
                                      const double& yaw) :
        weight_(weight),
        rpy_{roll, pitch, yaw}
    {
    }

private:
    const double weight_;
    const std::array<double, 3> rpy_;
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_3D_QUATERNION_HPP
