#ifndef CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
#define CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP

#include <cslibs_math_2d/linear/point.hpp>

#include <ceres/cost_function.h>
#include <ceres/cubic_interpolation.h>
#include <ceres/autodiff_cost_function.h>
#include <ceres/numeric_diff_cost_function.h>

namespace cslibs_ndt_2d {
namespace matching {

template <template <typename,typename> class child_t, typename base_t>
class CeresScanMatchCostFunctorCreator
{
public:
    template <typename points_t, typename ... args_t>
    static inline ceres::CostFunction* CreateAutoDiffCostFunction(double weight,
                                                                  points_t points,
                                                                  const args_t &...args)
    {
        using _child_t = child_t<base_t,points_t>;

        const auto count = points.size();
        return new ceres::AutoDiffCostFunction<_child_t, ceres::DYNAMIC, _child_t::N0, _child_t::N1>(
                    new _child_t(weight / std::sqrt(count), std::move(points), args...),
                    count);
    }

    template <ceres::NumericDiffMethodType method = ceres::FORWARD, typename points_t, typename ... args_t>
    static inline ceres::CostFunction* CreateNumericDiffCostFunction(double weight,
                                                                     points_t points,
                                                                     const args_t &...args)
    {
        using _child_t = child_t<base_t,points_t>;

        const auto count = points.size();
        return new ceres::NumericDiffCostFunction<_child_t, method, ceres::DYNAMIC, _child_t::N0, _child_t::N1>(
                    new _child_t(weight / std::sqrt(count), std::move(points), args...),
                    ceres::DO_NOT_TAKE_OWNERSHIP,
                    count);
    }
};

template <typename ndt_t>  // enable if with test for occupancy gridmap
class CeresDirectScanMatchCostFunctor
{
    using ivm_t = typename ndt_t::inverse_sensor_model_t;
    using point_t = typename ndt_t::point_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;

protected:
    explicit inline CeresDirectScanMatchCostFunctor(const ndt_t& map,
                                                    const typename ivm_t::Ptr& ivm) :
        map_(map),
        ivm_(ivm)
    {
    }

    inline void Evaluate(const double& x, const double& y, double* const value) const
    {
        *value = 1.0 - map_.sampleNonNormalized(point_t(x,y), ivm_);
    }

    template <typename JetT>
    inline void Evaluate(const JetT& x, const JetT& y, JetT* const value) const
    {
        const point_t pt(x.a, y.a);
        const bundle_t* bundle = map_.get(pt);

        const Eigen::Matrix<JetT, 2, 1> p(x,y);
        JetT retval(1);
        if (bundle) {
            for (std::size_t i=0; i<4; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;
                        auto sample = [&p,&di]() {
                            const auto &m = di->getMean();
                            const auto &i = di->getInformationMatrix();

                            const Eigen::Matrix<JetT, 2, 1> mean(JetT(m(0)), JetT(m(1)));
                            Eigen::Matrix<JetT, 2, 2> inf; inf << JetT(i(0,0)), JetT(i(0,1)), JetT(i(1,0)), JetT(i(1,1));
                            const Eigen::Matrix<JetT, 2, 1> q = p - mean;
                            const JetT exponent = -JetT(0.5) * q.transpose() * inf * q;
                            return ceres::exp(exponent);
                        };
                        retval -= 0.25 * bi->getOccupancy(ivm_) * sample();
                    }
                }
            }
        }
        *value = retval;
    }

private:
    const ndt_t& map_;
    const typename ivm_t::Ptr& ivm_;
};

template <typename ndt_t>  // enable if with test for occupancy gridmap
class CeresInterpolationScanMatchCostFunctor
{
    static constexpr int DATA_DIMENSION = 1;

    using ivm_t = typename ndt_t::inverse_sensor_model_t;
    using point_t = typename ndt_t::point_t;

    template <typename>
    friend class ceres::BiCubicInterpolator;

protected:
    explicit inline CeresInterpolationScanMatchCostFunctor(const ndt_t& map,
                                                           const typename ivm_t::Ptr& ivm,
                                                           const double& sampling_resolution) :
        map_(map),
        ivm_(ivm),
        sampling_resolution_(sampling_resolution),
        interpolator_(*this)
    {
    }

    template<typename T>
    inline void Evaluate(const T& x, const T& y, T* const value) const
    {
        interpolator_.Evaluate(x / sampling_resolution_,
                               y / sampling_resolution_,
                               value);
    }

private:
    inline void GetValue(const int row, const int column, double* const value) const
    {
        *value = 1.0 - map_.sampleNonNormalized(
                    point_t(row * sampling_resolution_,
                            column * sampling_resolution_),
                    ivm_);
    }

    const ndt_t& map_;
    const typename ivm_t::Ptr& ivm_;
    const double sampling_resolution_;
    const ceres::BiCubicInterpolator<CeresInterpolationScanMatchCostFunctor> interpolator_;
};
}
}

#endif // CSLIBS_NDT_2D_MATCHING_CERES_SCAN_MATCH_COST_FUNCTOR_HPP
