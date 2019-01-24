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
inline double validate(const double& val)
{
    if (val < 0.0 || val > 1.0 || !std::isnormal(val))
        return 0.0;
    return val;
}

inline void from(
        const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution)
{
    if (!src)
        return;
    src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::Gridmap;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const double bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [](const cslibs_math_2d::Point2d &p, const src_map_t::distribution_bundle_t &bundle) {
        return 0.25 * (validate(bundle.at(0)->data().sampleNonNormalized(p)) +
                       validate(bundle.at(1)->data().sampleNonNormalized(p)) +
                       validate(bundle.at(2)->data().sampleNonNormalized(p)) +
                       validate(bundle.at(3)->data().sampleNonNormalized(p)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2d p(bi[0] * bundle_resolution + k * sampling_resolution,
                                                bi[1] * bundle_resolution + l * sampling_resolution);
                dst->at((bi[0] - min_bi[0]) * chunk_step + k, (bi[1] - min_bi[1]) * chunk_step + l) = sample(p, b);
            }
        }
    });
}

inline void from(
        const cslibs_ndt_2d::static_maps::mono::Gridmap::Ptr  &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution)
{
    if (!src)
        return;

    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const std::size_t height = dst->getHeight();
    const std::size_t width  = dst->getWidth();
    const cslibs_ndt_2d::static_maps::mono::Gridmap::index_t min_index = src->getMinIndex();

    //// Todo: Check that code here!

    std::cerr << "converting " << width << " " << height << "\n";

    for(std::size_t i = 0 ; i < height ; ++i) {
        for(std::size_t j = 0 ; j < width ; ++j) {
            const cslibs_math_2d::Point2d  p(j * sampling_resolution,
                                             i * sampling_resolution);

            const cslibs_ndt_2d::static_maps::mono::Gridmap::index_t idx = {{min_index[0] + static_cast<int>(p(0) / src->getResolution()),
                                                                             min_index[1] + static_cast<int>(p(1) / src->getResolution())}};

            const double v = validate(src->sampleNonNormalized(src->getOrigin() * p,idx));
            if(v >= 0.0) {
                dst->at(j,i) = v;
            }
        }
    }
    std::cerr << "converted" << std::endl;
}

inline void from(
        const cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution,
        const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model)
{
    if (!src || !inverse_model)
        return;
    src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const double bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [&inverse_model](const cslibs_math_2d::Point2d &p, const src_map_t::distribution_bundle_t &bundle) {
        auto sample = [&p, &inverse_model](const src_map_t::distribution_t *d) {
            auto do_sample = [&p, &inverse_model, &d]() {
                return (d->getDistribution() && d->getDistribution()->valid()) ?
                            validate(d->getDistribution()->sampleNonNormalized(p)) *
                            validate(d->getOccupancy(inverse_model)) : 0.0;
            };
            return d ? do_sample() : 0.0;
        };
        return 0.25 * (sample(bundle.at(0)) +
                       sample(bundle.at(1)) +
                       sample(bundle.at(2)) +
                       sample(bundle.at(3)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2d p(bi[0] * bundle_resolution + k * sampling_resolution,
                                                bi[1] * bundle_resolution + l * sampling_resolution);
                dst->at((bi[0] - min_bi[0]) * chunk_step + k, (bi[1] - min_bi[1]) * chunk_step + l) = sample(p, b);
            }
        }
    });
}

inline void from(
        const cslibs_ndt_2d::dynamic_maps::WeightedOccupancyGridmap::Ptr &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution,
        const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model)
{
    if (!src || !inverse_model)
        return;
    src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_2d::dynamic_maps::WeightedOccupancyGridmap;
    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            std::ceil(src->getHeight() / sampling_resolution),
                            std::ceil(src->getWidth()  / sampling_resolution)));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const double bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto sample = [&inverse_model](const cslibs_math_2d::Point2d &p, const src_map_t::distribution_bundle_t &bundle) {
        auto sample = [&p, &inverse_model](const src_map_t::distribution_t *d) {
            auto do_sample = [&p, &inverse_model, &d]() {
                return (d->getDistribution() && d->getDistribution()->valid()) ?
                            validate(d->getDistribution()->sampleNonNormalized(p)) *
                            validate(d->getOccupancy(inverse_model)) : 0.0;
            };
            return d ? do_sample() : 0.0;
        };
        return 0.25 * (sample(bundle.at(0)) +
                       sample(bundle.at(1)) +
                       sample(bundle.at(2)) +
                       sample(bundle.at(3)));
    };

    using index_t = std::array<int, 2>;
    const index_t min_bi = src->getMinBundleIndex();

    src->traverse([&dst, &bundle_resolution, &sampling_resolution, &chunk_step, &min_bi, &sample]
                  (const index_t &bi, const src_map_t::distribution_bundle_t &b){
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                const cslibs_math_2d::Point2d p(bi[0] * bundle_resolution + k * sampling_resolution,
                                                bi[1] * bundle_resolution + l * sampling_resolution);
                dst->at((bi[0] - min_bi[0]) * chunk_step + k, (bi[1] - min_bi[1]) * chunk_step + l) = sample(p, b);
            }
        }
    });
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP
