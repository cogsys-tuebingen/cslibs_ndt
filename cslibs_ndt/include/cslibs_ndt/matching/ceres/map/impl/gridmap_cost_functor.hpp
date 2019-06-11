#ifndef CSLIBS_NDT_MATCHING_CERES_GRIDMAP_COST_FUNCTOR_HPP
#define CSLIBS_NDT_MATCHING_CERES_GRIDMAP_COST_FUNCTOR_HPP

#include <cslibs_ndt/map/map.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor.hpp>

#include <ceres/cubic_interpolation.h>
#include <type_traits>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class ScanMatchCostFunctor<
        cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>,
        Flag::DIRECT>
{
    using ndt_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>;

    using point_t = typename ndt_t::point_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;

protected:
    explicit inline ScanMatchCostFunctor(const ndt_t& map) :
        map_(map),
        resolution_inv_(1.0 / map_.getBundleResolution())
    {
        const auto& origin_inv = map_.getInitialOrigin().inverse();
        rot_ = Eigen::Rotation2D<double>(static_cast<double>(origin_inv.yaw())).toRotationMatrix();
        trans_ = Eigen::Matrix<double,2,1>(static_cast<double>(origin_inv.tx()),static_cast<double>(origin_inv.ty()));
    }

    template <int _D>
    inline void Evaluate(const Eigen::Matrix<double,_D,1>& q, double* const value) const
    {
        *value = 1.0 - map_.sampleNonNormalized(point_t(q(0),q(1)));
    }

    template <typename JetT, int _D>
    inline void Evaluate(const Eigen::Matrix<JetT,_D,1>& q, JetT* const value) const
    {
        const Eigen::Matrix<JetT,2,1>& p = q.template topRows<2>();

        const Eigen::Matrix<JetT,2,1>& p_prime = rot_ * p + trans_;
        const std::array<int,2> bi{{static_cast<int>(std::floor(p_prime(0).a * resolution_inv_)),
                                    static_cast<int>(std::floor(p_prime(1).a * resolution_inv_))}};

        const bundle_t* bundle = map_.get(bi);
        *value = JetT(1.0);
        if (bundle) {

            for (std::size_t i=0; i<ndt_t::bin_count; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const Eigen::Matrix<JetT,2,1> diff =
                                p_prime - di->getMean().template cast<double>();
                        const Eigen::Matrix<double,2,2> inf =
                                di->getInformationMatrix().template cast<double>();

                        const auto& sample = ::ceres::exp((-0.5 * diff.transpose() * inf * diff).value());
                        *value -= static_cast<double>(ndt_t::div_count) * sample;
                    }
                }
            }
        }
    }

private:
    const ndt_t& map_;

    const double resolution_inv_;
    Eigen::Matrix<double,2,2> rot_;
    Eigen::Matrix<double,2,1> trans_;
};

template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class ScanMatchCostFunctor<
        cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>,
        Flag::DIRECT>
{
    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>;

    using point_t = typename ndt_t::point_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;

protected:
    explicit inline ScanMatchCostFunctor(const ndt_t& map) :
        map_(map),
        resolution_inv_(1.0 / map_.getBundleResolution())
    {
        const auto& origin_inv = map_.getInitialOrigin().inverse();
        const auto& r = origin_inv.rotation();
        rot_ = Eigen::Quaternion<double>(
                    static_cast<double>(r.w()), static_cast<double>(r.x()),
                    static_cast<double>(r.y()), static_cast<double>(r.z()));
        trans_ = Eigen::Matrix<double,3,1>(
                    static_cast<double>(origin_inv.tx()), static_cast<double>(origin_inv.ty()), static_cast<double>(origin_inv.tz()));
    }

    template <int _D>
    inline void Evaluate(const Eigen::Matrix<double,_D,1>& q, double* const value) const
    {
        *value = 1.0 - map_.sampleNonNormalized(point_t(q(0),q(1),q(2)));
    }

    template <typename JetT, int _D>
    inline void Evaluate(const Eigen::Matrix<JetT,_D,1>& q, JetT* const value) const
    {
        const Eigen::Matrix<JetT,3,1>& p = q.template topRows<3>();

        const Eigen::Matrix<JetT,3,1>& p_prime = rot_ * p + trans_;
        const std::array<int,2> bi{{static_cast<int>(std::floor(p_prime(0).a * resolution_inv_)),
                                    static_cast<int>(std::floor(p_prime(1).a * resolution_inv_)),
                                    static_cast<int>(std::floor(p_prime(2).a * resolution_inv_))}};

        const bundle_t* bundle = map_.get(bi);
        *value = JetT(1.0);
        if (bundle) {

            for (std::size_t i=0; i<ndt_t::bin_count; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const Eigen::Matrix<JetT,3,1> diff =
                                p_prime - di->getMean().template cast<double>();
                        const Eigen::Matrix<double,3,3> inf =
                                di->getInformationMatrix().template cast<double>();

                        const auto& sample = ::ceres::exp((-0.5 * diff.transpose() * inf * diff).value());
                        *value -= static_cast<double>(ndt_t::div_count) * sample;
                    }
                }
            }
        }
    }

private:
    const ndt_t& map_;

    const double resolution_inv_;
    Eigen::Quaternion<double> rot_;
    Eigen::Matrix<double,3,1> trans_;
};

// only possible for maps of dimension 2
template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class ScanMatchCostFunctor<
        cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>,
        Flag::INTERPOLATION>
{
    using ndt_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>;

    using point_t = typename ndt_t::point_t;

    using transform_t = typename ndt_t::pose_t;
    using bundle_t = typename ndt_t::distribution_bundle_t;
    using index_t = typename ndt_t::index_t;

    template <typename>
    friend class ::ceres::BiCubicInterpolator;

    static constexpr int DATA_DIMENSION = 1;

protected:
    explicit inline ScanMatchCostFunctor(const ndt_t& map,
                                         const double& sampling_resolution) :
        map_(map),
        sampling_resolution_(sampling_resolution),
        interpolator_(*this)
    {
    }

    template <int _D>
    inline void Evaluate(const Eigen::Matrix<double,_D,1>& q, double* const value) const
    {
        *value = 1.0 - map_.sampleNonNormalized(point_t(q(0),q(1)));
    }

    template <typename JetT, int _D>
    inline void Evaluate(const Eigen::Matrix<JetT,_D,1>& q, JetT* const value) const
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
                            column * sampling_resolution_));
    }

    const ndt_t& map_;
    const double sampling_resolution_;
    const ::ceres::BiCubicInterpolator<ScanMatchCostFunctor<ndt_t,Flag::INTERPOLATION>> interpolator_;
};

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_GRIDMAP_COST_FUNCTOR_HPP
