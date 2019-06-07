#ifndef CSLIBS_NDT_MATCHING_CERES_OCCUPANCY_GRIDMAP_COST_FUNCTOR_HPP
#define CSLIBS_NDT_MATCHING_CERES_OCCUPANCY_GRIDMAP_COST_FUNCTOR_HPP

#include <cslibs_ndt/map/map.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor.hpp>

#include <ceres/cubic_interpolation.h>
#include <type_traits>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <cslibs_ndt::map::tags::option option_t,
          std::size_t Dim,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class ScanMatchCostFunctor<
        cslibs_ndt::map::Map<option_t,Dim,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>,
        Flag::DIRECT>
{
    using ndt_t = cslibs_ndt::map::Map<option_t,Dim,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>;

    using ivm_t = typename ndt_t::inverse_sensor_model_t;
    using point_t = typename ndt_t::point_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;

protected:
    explicit inline ScanMatchCostFunctor(const ndt_t& map,
                                         const typename ivm_t::Ptr& ivm) :
        map_(map),
        ivm_(ivm)
    {
    }

    template <typename JetT>
    inline Eigen::Matrix<JetT,2,1> transformToMap(const Eigen::Matrix<JetT,2,1>& p) const
    {
        static const auto& origin_inv = map_.getInitialOrigin().inverse();

        static const Eigen::Matrix<JetT,2,2> rot =
                Eigen::Rotation2D<JetT>(JetT(origin_inv.yaw())).toRotationMatrix();
        static const Eigen::Matrix<JetT,2,1> trans(JetT(origin_inv.tx()),JetT(origin_inv.ty()));

        return rot * p + trans;
    }

    template <typename JetT>
    inline Eigen::Matrix<JetT,3,1> transformToMap(const Eigen::Matrix<JetT,3,1>& p) const
    {
        static const auto& origin_inv = map_.getInitialOrigin().inverse();

        static const auto& r = origin_inv.rotation();
        static const Eigen::Quaternion<JetT> rot(JetT(r.w()), JetT(r.x()), JetT(r.y()), JetT(r.z()));
        static const Eigen::Matrix<JetT,3,1> trans(JetT(origin_inv.tx()), JetT(origin_inv.ty()), JetT(origin_inv.tz()));

        return rot * p + trans;
    }

    template <int _D>
    inline void Evaluate(const Eigen::Matrix<double,_D,1>& q, double* const value) const
    {
        Eigen::Matrix<_T,Dim,1> p = Eigen::Matrix<_T,Dim,1>::Zero();
        for (std::size_t i=0; i<std::min(_D,static_cast<int>(Dim)); ++i)
            p(i) = q(i);

        *value = 1.0 - static_cast<double>(map_.sampleNonNormalized(point_t(p), ivm_));
    }

    template <typename JetT, int _D>
    inline void Evaluate(const Eigen::Matrix<JetT,_D,1>& q, JetT* const value) const
    {
        point_t pt;
        Eigen::Matrix<JetT,Dim,1> p = Eigen::Matrix<JetT,Dim,1>::Zero();
        for (std::size_t i=0; i<std::min(_D,static_cast<int>(Dim)); ++i) {
            p(i)  = q(i);
            pt(i) = q(i).a;
        }

        const bundle_t* bundle = map_.get(pt);
        *value = static_cast<JetT>(1.0);
        if (bundle) {
            const Eigen::Matrix<JetT,Dim,1>& p_prime = transformToMap(p);

            for (std::size_t i=0; i<ndt_t::bin_count; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    const double occ = static_cast<double>(bi->getOccupancy(ivm_));
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const Eigen::Matrix<JetT,Dim,1> diff =
                                p_prime - di->getMean().template cast<double>();
                        const Eigen::Matrix<double,Dim,Dim> inf =
                                di->getInformationMatrix().template cast<double>();

                        const JetT sample = ::ceres::exp((-0.5 * diff.transpose() * inf * diff).eval().value());
                        *value -= static_cast<double>(ndt_t::div_count) * sample * occ;
                    }
                }
            }
        }
    }

private:
    const ndt_t& map_;
    const typename ivm_t::Ptr& ivm_;
};

// only possible for maps of dimension 2
template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class ScanMatchCostFunctor<
        cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>,
        Flag::INTERPOLATION>
{
    using ndt_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>;

    using ivm_t = typename ndt_t::inverse_sensor_model_t;
    using point_t = typename ndt_t::point_t;

    using transform_t = typename ndt_t::pose_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;
    using index_t = typename ndt_t::index_t;

    template <typename>
    friend class ::ceres::BiCubicInterpolator;

    static constexpr int DATA_DIMENSION = 1;

protected:
    explicit inline ScanMatchCostFunctor(const ndt_t& map,
                                         const typename ivm_t::Ptr& ivm,
                                         const double& sampling_resolution) :
        map_(map),
        ivm_(ivm),
        sampling_resolution_(sampling_resolution),
        interpolator_(*this)
    {
    }

    template <typename T, int _D>
    inline void Evaluate(const Eigen::Matrix<T,_D,1>& q, T* const value) const
    {
        interpolator_.Evaluate(q(0) / sampling_resolution_,
                               q(1) / sampling_resolution_,
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
    const ::ceres::BiCubicInterpolator<ScanMatchCostFunctor<ndt_t,Flag::INTERPOLATION>> interpolator_;
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_OCCUPANCY_GRIDMAP_COST_FUNCTOR_HPP
