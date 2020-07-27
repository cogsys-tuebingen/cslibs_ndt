#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP

#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_creator.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor.hpp>

#include <cslibs_math_3d/linear/point.hpp>
#include <cslibs_math_3d/linear/quaternion.hpp>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <typename base_t, typename points_t>
class ScanMatchCostFunctor3dRPY : public base_t
{
    static constexpr int N0 = 3;
    static constexpr int N1 = 3;

    template <template <typename,typename> class, typename>
    friend class ceres::ScanMatchCostFunctorCreator;

public:
    template <typename T>
    inline Eigen::Quaternion<T> toEigen(const T* const raw_rotation_rpy) const
    {
        const T roll_2    = raw_rotation_rpy[0] * 0.5;
        const T pitch_2   = raw_rotation_rpy[1] * 0.5;
        const T yaw_2     = raw_rotation_rpy[2] * 0.5;
        const T cos_roll  = ::ceres::cos(roll_2);
        const T sin_roll  = ::ceres::sin(roll_2);
        const T cos_pitch = ::ceres::cos(pitch_2);
        const T sin_pitch = ::ceres::sin(pitch_2);
        const T cos_yaw   = ::ceres::cos(yaw_2);
        const T sin_yaw   = ::ceres::sin(yaw_2);

        const T data_x = sin_roll * cos_pitch * cos_yaw - cos_roll * sin_pitch * sin_yaw;
        const T data_y = cos_roll * sin_pitch * cos_yaw + sin_roll * cos_pitch * sin_yaw;
        const T data_z = cos_roll * cos_pitch * sin_yaw - sin_roll * sin_pitch * cos_yaw;
        const T data_w = cos_roll * cos_pitch * cos_yaw + sin_roll * sin_pitch * sin_yaw;

        return Eigen::Quaternion<T>{ data_w, data_x, data_y, data_z };
    }

    template <typename T>
    inline bool operator()(const T* const raw_translation, const T* const raw_rotation_rpy, T* residual) const
    {
        const Eigen::Matrix<T, 3, 1> translation(raw_translation[0], raw_translation[1], raw_translation[2]);
        const Eigen::Quaternion<T> rotation = toEigen(raw_rotation_rpy);

        std::size_t i = 0;
        //const double size = static_cast<double>(points_.size());
        for (const auto& point : points_) {
            const Eigen::Matrix<T, 3, 1> local(T(point(0)), T(point(1)), T(point(2)));
            const Eigen::Matrix<T, 3, 1> in_world = rotation * local + translation;

            this->Evaluate(in_world, &residual[i]);
            residual[i] = weight_ * residual[i]; //::ceres::sqrt(weight_) * residual[i] / size;
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

template <typename ndt_t, Flag flag_t>
using ScanMatchCostFunctor3dRPYCreator =
ScanMatchCostFunctorCreator<ScanMatchCostFunctor3dRPY, ScanMatchCostFunctor<ndt_t, flag_t>>;

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_3D_RPY_HPP
