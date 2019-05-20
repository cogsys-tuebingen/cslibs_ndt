#ifndef CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP
#define CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP

#include <cslibs_ndt_2d/matching/ceres/scan_match_cost_functor.hpp>
#include <cslibs_math_3d/linear/point.hpp>
#include <cslibs_math_3d/linear/quaternion.hpp>

namespace cslibs_ndt_2d {
namespace matching {
namespace ceres {

template <typename base_t, typename points_t>
class ScanMatchCostFunctor3dRPY : public base_t
{
    static constexpr int N0 = 3;
    static constexpr int N1 = 3;

    template <template <typename,typename> class, typename>
    friend class ScanMatchCostFunctorCreator;

public:
    template<typename T>
    inline bool operator()(const T* const raw_translation, const T* const raw_rotation_rpy, T* residual) const
    {
        const Eigen::Matrix<T, 3, 1> translation(raw_translation[0], raw_translation[1], raw_translation[2]);
        const Eigen::Quaternion<T>   rotation =
                cslibs_math_3d::Quaternion<T>(T(raw_rotation_rpy[0]),
                                              T(raw_rotation_rpy[1]),
                                              T(raw_rotation_rpy[2])).toEigen();

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
    explicit inline ScanMatchCostFunctor3dRPY(const double& weight,
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
using DirectScanMatchCostFunctor3dRPY =
ScanMatchCostFunctorCreator<ScanMatchCostFunctor3dRPY, DirectScanMatchCostFunctor<ndt_t>>;

template <typename ndt_t>
using InterpolationScanMatchCostFunctor3dRPY =
ScanMatchCostFunctorCreator<ScanMatchCostFunctor3dRPY, InterpolationScanMatchCostFunctor<ndt_t>>;

}
}
}

#endif // CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP
