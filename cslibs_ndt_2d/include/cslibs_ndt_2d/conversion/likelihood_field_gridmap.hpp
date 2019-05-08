#ifndef CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP
#define CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP

#include <cslibs_ndt/map/map.hpp>

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_gridmaps/static_maps/likelihood_field_gridmap.h>
#include <cslibs_gridmaps/static_maps/algorithms/distance_transform.hpp>

namespace cslibs_ndt_2d {
namespace conversion {
template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
inline void from(
        const typename cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,T,backend_t,dynamic_backend_t> &src,
        typename cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T,T>::Ptr &dst,
        const T &sampling_resolution,
        const T &maximum_distance = 2.0,
        const T &sigma_hit        = 0.5,
        const T &threshold        = 0.169,
        const bool allocate_all   = true)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    assert(threshold <= 1.0);
    assert(threshold >= 0.0);
    const T exp_factor_hit = (cslibs_math::utility::traits<T>::Half / (sigma_hit * sigma_hit));

    using src_map_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,T,backend_t,dynamic_backend_t>;
    using dst_map_t = cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T,T>;
    dst.reset(new dst_map_t(src.getOrigin(),
                            sampling_resolution,
                            std::ceil(src.getHeight() / sampling_resolution),
                            std::ceil(src.getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0.5);//T());

    const T bundle_resolution = src.getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [](const cslibs_math_2d::Point2<T> &p, const typename src_map_t::distribution_bundle_t &bundle) {
        return src_map_t::div_count * (bundle.at(0)->data().sampleNonNormalized(p) +
                                       bundle.at(1)->data().sampleNonNormalized(p) +
                                       bundle.at(2)->data().sampleNonNormalized(p) +
                                       bundle.at(3)->data().sampleNonNormalized(p));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src.getMinBundleIndex();

    const auto& origin = src.getInitialOrigin();
    src.traverse([&dst, &origin, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const typename src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2<T> p(static_cast<T>(bi[0]) * bundle_resolution + static_cast<T>(k) * sampling_resolution,
                                                  static_cast<T>(bi[1]) * bundle_resolution + static_cast<T>(l) * sampling_resolution);
                const std::size_t u = (bi[0] - min_bi[0]) * chunk_step + k;
                const std::size_t v = (bi[1] - min_bi[1]) * chunk_step + l;
                dst->at(u,v) = sample(origin * p, b);
            }
        }
    });

    std::vector<T> occ = dst->getData();
    cslibs_gridmaps::static_maps::algorithms::DistanceTransform<T,T,T> distance_transform(
                sampling_resolution, maximum_distance, threshold);
    distance_transform.apply(occ, dst->getWidth(), dst->getData());

    std::for_each(dst->getData().begin(),
                  dst->getData().end(),
                  [&exp_factor_hit] (T &z) {z = std::exp(-z * z * exp_factor_hit);});
}

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
inline void from(
        const typename cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,T,backend_t,dynamic_backend_t> &src,
        typename cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T,T>::Ptr &dst,
        const T &sampling_resolution,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &inverse_model,
        const T &maximum_distance = 2.0,
        const T &sigma_hit        = 0.5,
        const T &threshold        = 0.169,
        const bool allocate_all   = true)
{
    if (!inverse_model)
        return;
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    assert(threshold <= 1.0);
    assert(threshold >= 0.0);
    const T exp_factor_hit = (cslibs_math::utility::traits<T>::Half / (sigma_hit * sigma_hit));

    using src_map_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,T,backend_t,dynamic_backend_t>;
    using dst_map_t = cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T>;
    dst.reset(new dst_map_t(src.getOrigin(),
                            sampling_resolution,
                            std::ceil(src.getHeight() / sampling_resolution),
                            std::ceil(src.getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0.5);//T());

    const T bundle_resolution = src.getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [&inverse_model](const cslibs_math_2d::Point2<T> &p, const typename src_map_t::distribution_bundle_t &bundle) {
        auto sample = [&p, &inverse_model](const typename src_map_t::distribution_t *d) {
            auto do_sample = [&p, &inverse_model, &d]() {
                const auto &handle = d;
                return handle->getDistribution() ?
                            handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(inverse_model) : T();
            };
            return d ? do_sample() : T();
        };
        return src_map_t::div_count * (sample(bundle.at(0)) +
                                       sample(bundle.at(1)) +
                                       sample(bundle.at(2)) +
                                       sample(bundle.at(3)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src.getMinBundleIndex();

    const auto& origin = src.getInitialOrigin();
    src.traverse([&dst, &origin, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const typename src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2<T> p(static_cast<T>(bi[0]) * bundle_resolution + static_cast<T>(k) * sampling_resolution,
                                                  static_cast<T>(bi[1]) * bundle_resolution + static_cast<T>(l) * sampling_resolution);
                const std::size_t u = (bi[0] - min_bi[0]) * chunk_step + k;
                const std::size_t v = (bi[1] - min_bi[1]) * chunk_step + l;
                dst->at(u,v) = sample(origin * p, b);
            }
        }
    });

    std::vector<T> occ = dst->getData();
    cslibs_gridmaps::static_maps::algorithms::DistanceTransform<T,T,T> distance_transform(
                sampling_resolution, maximum_distance, threshold);
    distance_transform.apply(occ, dst->getWidth(), dst->getData());

    std::for_each(dst->getData().begin(),
                  dst->getData().end(),
                  [&exp_factor_hit] (T &z) {z = std::exp(-z * z * exp_factor_hit);});
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::dynamic_maps::Gridmap<T>::Ptr &src,
        typename cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T,T>::Ptr &dst,
        const T &sampling_resolution,
        const T &maximum_distance = 2.0,
        const T &sigma_hit        = 0.5,
        const T &threshold        = 0.169,
        const bool allocate_all   = true)
{
    if (!src)
        return;
    return from<
            cslibs_ndt::map::tags::dynamic_map,
            T,
            cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::default_backend_t,
            cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::default_dynamic_backend_t>(
                *src, dst, sampling_resolution, maximum_distance, sigma_hit, threshold, allocate_all);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        typename cslibs_gridmaps::static_maps::LikelihoodFieldGridmap<T,T>::Ptr &dst,
        const T &sampling_resolution,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &inverse_model,
        const T &maximum_distance = 2.0,
        const T &sigma_hit        = 0.5,
        const T &threshold        = 0.169,
        const bool allocate_all   = true)
{
    if (!src)
        return;
    return from<
            cslibs_ndt::map::tags::dynamic_map,
            T,
            cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::default_backend_t,
            cslibs_ndt::map::tags::default_types<cslibs_ndt::map::tags::dynamic_map>::default_dynamic_backend_t>(
                *src, dst, sampling_resolution, inverse_model, maximum_distance, sigma_hit, threshold, allocate_all);
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP
