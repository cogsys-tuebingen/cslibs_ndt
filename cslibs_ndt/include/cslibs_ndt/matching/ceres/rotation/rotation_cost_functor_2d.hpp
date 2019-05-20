#ifndef CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_2D_HPP
#define CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_2D_HPP

#include <ceres/cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/jet.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

class RotationCostFunctor2d
{
public:
    static ::ceres::CostFunction* CreateAutoDiffCostFunction(double weight,
                                                             const double& rotation)
    {
        return new ::ceres::AutoDiffCostFunction<RotationCostFunctor2d, 1, 1>(
                new RotationCostFunctor2d(weight, rotation)
        );
    }

    template<typename T>
    bool operator()(const T* const rotation, T* residual) const
    {
        T delta;
        T r0(rotation_);
        AngleDifference(&r0, &(rotation[0]), &delta);

        residual[0] = weight_ * delta;
        return true;
    }

private:
    template<typename T>
    static void AngleDifference(const T* const r0, const T* const r, T* const d)
    {
        static const double _2_M_PI = 2.0 * M_PI;
        auto norm = [](const T* const r) {
            return ::ceres::atan2(::ceres::sin(*r), ::ceres::cos(*r));
        };

        const auto a = norm(r0);
        const auto b = norm(r);

        const auto d1 = a - b;
        const auto d2 = (_2_M_PI - ::ceres::abs(d1)) * (d1 > T(0) ? T(-1) : T(1));
        *d = ::ceres::abs(d1) < ::ceres::abs(d2) ? d1 : d2;
    }

    explicit RotationCostFunctor2d(double weight,
                                   const double& rotation) :
        weight_(weight),
        rotation_(rotation)
    {
    }

private:
    const double weight_;
    const double rotation_;
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_ROTATION_COST_FUNCTOR_2D_HPP
