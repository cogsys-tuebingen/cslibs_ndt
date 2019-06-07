#ifndef CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_2D_HPP
#define CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_2D_HPP

#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_creator.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor.hpp>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <typename base_t, typename points_t>
class ScanMatchCostFunctor2d : public base_t
{
    static constexpr int N0 = 2;
    static constexpr int N1 = 1;

    template <template <typename,typename> class, typename>
    friend class ScanMatchCostFunctorCreator;

public:
    template<typename T>
    inline bool operator()(const T* const raw_translation, const T* const raw_rotation, T* residual) const
    {
        const Eigen::Matrix<T,2,1> translation(raw_translation[0], raw_translation[1]);
        const Eigen::Matrix<T,2,2> rotation =
                Eigen::Rotation2D<T>(raw_rotation[0]).toRotationMatrix();

        std::size_t i = 0;
        const double size = static_cast<double>(points_.size());
        for (const auto& point : points_) {
            const Eigen::Matrix<T, 2, 1> local(T(point(0)), T(point(1)));
            const Eigen::Matrix<T, 2, 1> in_world = rotation * local + translation;

            this->Evaluate(in_world, &residual[i]);
            residual[i] = ::ceres::sqrt(weight_) * residual[i] / size;
            ++i;
        }
        return true;
    }

private:
    template <typename ... args_t>
    explicit inline ScanMatchCostFunctor2d(const double& weight,
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
using ScanMatchCostFunctor2dCreator =
ScanMatchCostFunctorCreator<ScanMatchCostFunctor2d, ScanMatchCostFunctor<ndt_t, flag_t>>;

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_2D_HPP
