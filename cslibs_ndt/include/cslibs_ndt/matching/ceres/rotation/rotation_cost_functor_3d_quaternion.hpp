#ifndef CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_QUATERNION_3D_HPP
#define CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_QUATERNION_3D_HPP

#include <cslibs_math_3d/linear/quaternion.hpp>

#include <ceres/cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/jet.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

class RotationCostFunctor3dQuaternion
{
public:
    static ::ceres::CostFunction* CreateAutoDiffCostFunction(const double& weight,
                                                             const cslibs_math_3d::Quaterniond& rotation)
    {
        return new ::ceres::AutoDiffCostFunction<RotationCostFunctor3dQuaternion, 4, 4>(
                new RotationCostFunctor3dQuaternion(weight, rotation)
        );
    }

    template<typename T>
    bool operator()(const T* const rotation_wxyz, T* residual) const
    {
        std::array<T, 4> delta;
        QuaternionProduct(rotation_inverse_wxyz_.data(), rotation_wxyz, delta.data());

        residual[0] = ::ceres::sqrt(weight_) * ::ceres::abs(delta[0]);
        residual[1] = ::ceres::sqrt(weight_) * ::ceres::abs(delta[1]);
        residual[2] = ::ceres::sqrt(weight_) * ::ceres::abs(delta[2]);
        residual[3] = ::ceres::sqrt(weight_) * ::ceres::abs(delta[3]);
        return true;
    }

private:
    template<typename T>
    static void QuaternionProduct(const double* const z, const T* const w, T* const zw)
    {
        zw[0] = z[0] * w[0] - z[1] * w[1] - z[2] * w[2] - z[3] * w[3];
        zw[1] = z[0] * w[1] + z[1] * w[0] + z[2] * w[3] - z[3] * w[2];
        zw[2] = z[0] * w[2] - z[1] * w[3] + z[2] * w[0] + z[3] * w[1];
        zw[3] = z[0] * w[3] + z[1] * w[2] - z[2] * w[1] + z[3] * w[0];
    }

    explicit RotationCostFunctor3dQuaternion(const double& weight,
                                             const cslibs_math_3d::Quaterniond& rotation) :
        weight_(weight),
        rotation_inverse_wxyz_({
            rotation.w(),
            -rotation.x(),
            -rotation.y(),
            -rotation.z()})
    {
    }

private:
    const double weight_;
    const std::array<double, 4> rotation_inverse_wxyz_;
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_3D_QUATERNION_HPP
