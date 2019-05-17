#ifndef CSLIBS_NDT_MATCHING_CERES_TRANSLATION_COST_FUNCTOR_2D_HPP
#define CSLIBS_NDT_MATCHING_CERES_TRANSLATION_COST_FUNCTOR_2D_HPP

#include <cslibs_math_2d/linear/vector.hpp>

#include <ceres/cost_function.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/jet.h>

namespace cslibs_ndt {
namespace matching {

class CeresTranslationCostFunctor2d
{
public:
    static ceres::CostFunction* CreateAutoDiffCostFunction(double weight,
                                                           const cslibs_math_2d::Vector2d& translation)
    {
        return new ceres::AutoDiffCostFunction<CeresTranslationCostFunctor2d, 2, 2>(
                new CeresTranslationCostFunctor2d(weight, translation)
        );
    }

    template<typename T>
    bool operator()(const T* const translation, T* residual) const
    {
        for (auto i = 0; i < 2; ++i)
            residual[i] = weight_ * (translation[i] - translation_[i]);

        return true;
    }

private:
    explicit CeresTranslationCostFunctor2d(double weight,
                                           const cslibs_math_2d::Vector2d& translation) :
        weight_(weight),
        translation_({
                translation(0),
                translation(1)})
    {
    }

private:
    const double weight_;
    const std::array<double, 2> translation_;
};

}
}

#endif // CSLIBS_NDT_MATCHING_CERES_TRANSLATION_COST_FUNCTOR_2D_HPP
