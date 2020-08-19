#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_QUATERNION_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_QUATERNION_HPP

#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_creator.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor.hpp>

#include <cslibs_math_3d/linear/point.hpp>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <typename base_t, typename points_t>
class ScanMatchCostFunctor3dQuaternion : public base_t
{
    static constexpr int N0 = 3;
    static constexpr int N1 = 4;

    template <template <typename,typename> class, typename>
    friend class ScanMatchCostFunctorCreator;

public:
    template<typename T>
    inline bool operator()(const T* const raw_translation, const T* const raw_rotation_wxyz, T* residual) const
    {
        const Eigen::Matrix<T, 3, 1> translation(raw_translation[0], raw_translation[1], raw_translation[2]);
        const Eigen::Quaternion<T>   rotation(raw_rotation_wxyz[0], raw_rotation_wxyz[1], raw_rotation_wxyz[2], raw_rotation_wxyz[3]);

        std::size_t i = 0;
        //const double size = static_cast<double>(points_.size());
        for (const auto& point : points_) {
            const Eigen::Matrix<T, 3, 1> local(T(point(0)), T(point(1)), T(point(2)));
            const Eigen::Matrix<T, 3, 1> in_world = rotation * local + translation;

            this->Evaluate(in_world, &residual[i]);            
            if (residual[i] == -residual[i]) // only nan test that works
                residual[i] = T(0.);
            residual[i] = weight_ * residual[i]; //::ceres::sqrt(weight_) * residual[i] / size;
            ++i;
        }
        return true;
    }

private:
    template <typename ... args_t>
    explicit inline ScanMatchCostFunctor3dQuaternion(const double& weight,
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

template <typename ndt_t, Flag flag_t>
using ScanMatchCostFunctor3dQuaternionCreator =
ScanMatchCostFunctorCreator<ScanMatchCostFunctor3dQuaternion, ScanMatchCostFunctor<ndt_t, flag_t>>;

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_QUATERNION_HPP
