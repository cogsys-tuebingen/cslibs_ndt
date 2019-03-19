#ifndef CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/weighted_occupancy_gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/mono_gridmap.hpp>

#include <cslibs_ndt_2d/conversion/gridmap.hpp>
#include <cslibs_ndt_2d/conversion/occupancy_gridmap.hpp>

#include <cslibs_gridmaps/static_maps/probability_gridmap.h>

namespace cslibs_ndt_2d {
namespace conversion {
template <typename T>
inline T validate(const T& val)
{
    if (val < 0 || val > cslibs_math::utility::traits<T>::One || !std::isnormal(val))
        return T();
    return val;
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::dynamic_maps::Gridmap<T>::Ptr &src,
        typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr &dst,
        const T sampling_resolution,
        const bool allocate_all = true)
{
    if (!src)
        return;
    if (allocate_all)
        src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::Gridmap<T>;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), T());

    const T bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [](const cslibs_math_2d::Point2<T> &p, const typename src_map_t::distribution_bundle_t &bundle) {
        return src_map_t::div_count * (validate(bundle.at(0)->data().sampleNonNormalized(p)) +
                                       validate(bundle.at(1)->data().sampleNonNormalized(p)) +
                                       validate(bundle.at(2)->data().sampleNonNormalized(p)) +
                                       validate(bundle.at(3)->data().sampleNonNormalized(p)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&src, &dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const typename src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2<T> p(static_cast<T>(bi[0]) * bundle_resolution + static_cast<T>(k) * sampling_resolution,
                                                  static_cast<T>(bi[1]) * bundle_resolution + static_cast<T>(l) * sampling_resolution);
                const std::size_t u = (bi[0] - min_bi[0]) * chunk_step + k;
                const std::size_t v = (bi[1] - min_bi[1]) * chunk_step + l;
                dst->at(u,v) = sample(p, b);
            }
        }
    });
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::static_maps::mono::Gridmap<T>::Ptr  &src,
        typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr &dst,
        const T sampling_resolution)
{
    if (!src)
        return;

    using src_map_t = cslibs_ndt_2d::static_maps::mono::Gridmap<T>;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), T());

    const std::size_t height = dst->getHeight();
    const std::size_t width  = dst->getWidth();
    const typename src_map_t::index_t min_index = src->getMinIndex();

    //// Todo: Check that code here!

    std::cerr << "converting " << width << " " << height << "\n";

    for(std::size_t i = 0 ; i < height ; ++i) {
        for(std::size_t j = 0 ; j < width ; ++j) {
            const cslibs_math_2d::Point2<T> p(j * sampling_resolution,
                                              i * sampling_resolution);

            const typename src_map_t::index_t idx = {{min_index[0] + static_cast<int>(p(0) / src->getResolution()),
                                                      min_index[1] + static_cast<int>(p(1) / src->getResolution())}};

            const T v = validate(src->sampleNonNormalized(src->getOrigin() * p,idx));
            if (v >= 0.0) {
                dst->at(j,i) = v;
            }
        }
    }
    std::cerr << "converted" << std::endl;
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr &dst,
        const T sampling_resolution,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &inverse_model,
        const bool allocate_all = true)
{
    if (!src || !inverse_model)
        return;
    if (allocate_all)
        src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), T());

    const T bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [&inverse_model](const cslibs_math_2d::Point2<T> &p, const typename src_map_t::distribution_bundle_t &bundle) {
        auto sample = [&p, &inverse_model](const typename src_map_t::distribution_t *d) {
            auto do_sample = [&p, &inverse_model, &d]() {
                return (d->getDistribution() && d->getDistribution()->valid()) ?
                            validate(d->getDistribution()->sampleNonNormalized(p)) *
                            validate(d->getOccupancy(inverse_model)) : T();
            };
            return d ? do_sample() : T();
        };
        return src_map_t::div_count * (sample(bundle.at(0)) +
                                       sample(bundle.at(1)) +
                                       sample(bundle.at(2)) +
                                       sample(bundle.at(3)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const typename src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2<T> p(static_cast<T>(bi[0]) * bundle_resolution + static_cast<T>(k) * sampling_resolution,
                                                  static_cast<T>(bi[1]) * bundle_resolution + static_cast<T>(l) * sampling_resolution);
                const std::size_t u = (bi[0] - min_bi[0]) * chunk_step + k;
                const std::size_t v = (bi[1] - min_bi[1]) * chunk_step + l;
                dst->at(u,v) = sample(p, b);
            }
        }
    });
}

template <typename T>
inline void from(
        const typename cslibs_ndt_2d::dynamic_maps::WeightedOccupancyGridmap<T>::Ptr &src,
        typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr &dst,
        const T sampling_resolution,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &inverse_model,
        const bool allocate_all = true)
{
    if (!src || !inverse_model)
        return;
    if (allocate_all)
        src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::WeightedOccupancyGridmap<T>;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), T());

    const T bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [&inverse_model](const cslibs_math_2d::Point2<T> &p, const typename src_map_t::distribution_bundle_t &bundle) {
        auto sample = [&p, &inverse_model](const typename src_map_t::distribution_t *d) {
            auto do_sample = [&p, &inverse_model, &d]() {
                return (d->getDistribution() && d->getDistribution()->valid()) ?
                            validate(d->getDistribution()->sampleNonNormalized(p)) *
                            validate(d->getOccupancy(inverse_model)) : T();
            };
            return d ? do_sample() : T();
        };
        return src_map_t::div_count * (sample(bundle.at(0)) +
                                       sample(bundle.at(1)) +
                                       sample(bundle.at(2)) +
                                       sample(bundle.at(3)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const typename src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2<T> p(static_cast<T>(bi[0]) * bundle_resolution + static_cast<T>(k) * sampling_resolution,
                                                  static_cast<T>(bi[1]) * bundle_resolution + static_cast<T>(l) * sampling_resolution);
                const std::size_t u = (bi[0] - min_bi[0]) * chunk_step + k;
                const std::size_t v = (bi[1] - min_bi[1]) * chunk_step + l;
                dst->at(u,v) = sample(p, b);
            }
        }
    });
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP
