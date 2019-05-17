#ifndef CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_HPP
#define CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_HPP

#include <cslibs_ndt_2d/matching/ceres_scan_match_cost_functor.hpp>
#include <cslibs_math_3d/linear/point.hpp>

namespace cslibs_ndt_2d {
namespace matching {

template <typename base_t, typename points_t>
class CeresScanMatchCostFunctor3d : public base_t
{
    static constexpr int N0 = 3;
    static constexpr int N1 = 4;

    template <template <typename,typename> class, typename>
    friend class CeresScanMatchCostFunctorCreator;

public:
    template<typename T>
    inline bool operator()(const T* const raw_translation, const T* const raw_rotation_wxyz, T* residual) const
    {
        const Eigen::Matrix<T, 3, 1> translation(raw_translation[0], raw_translation[1], raw_translation[2]);
        const Eigen::Quaternion<T>   rotation(raw_rotation_wxyz[0], raw_rotation_wxyz[1], raw_rotation_wxyz[2], raw_rotation_wxyz[3]);

        std::size_t i = 0;
        for (const auto& point : points_) {
            const Eigen::Matrix<T, 3, 1> local(T(point(0)), T(point(1)), T(point(2)));
            const Eigen::Matrix<T, 3, 1> in_world = rotation * local + translation;

            this->Evaluate(in_world(0), in_world(1), &residual[i]);
            residual[i] = weight_ * residual[i];
            ++i;
        }
        return true;
    }

private:
    template <typename ... args_t>
    explicit inline CeresScanMatchCostFunctor3d(const double& weight,
                                                points_t&& points,
                                                const args_t &...args) :
        base_t(args...),
        weight_(weight),
        points_(points)
    {
    }

    const double weight_;
    const points_t points_;
};

template <typename ndt_t>
using CeresDirectScanMatchCostFunctor3d =
CeresScanMatchCostFunctorCreator<CeresScanMatchCostFunctor3d, CeresDirectScanMatchCostFunctor<ndt_t>>;

template <typename ndt_t>
using CeresInterpolationScanMatchCostFunctor3d =
CeresScanMatchCostFunctorCreator<CeresScanMatchCostFunctor3d, CeresInterpolationScanMatchCostFunctor<ndt_t>>;

}
}

#endif // CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_HPP
