#ifndef CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/flat_gridmap.hpp>

#include <cslibs_ndt_2d/conversion/gridmap.hpp>
#include <cslibs_ndt_2d/conversion/occupancy_gridmap.hpp>

#include <cslibs_gridmaps/static_maps/probability_gridmap.h>

namespace cslibs_ndt_2d {
namespace conversion {
inline void from(
        const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution)
{
    if (!src)
        return;

    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            src->getHeight() / sampling_resolution,
                            src->getWidth()  / sampling_resolution));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const double bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinBundleIndex();
    const index_t max_distribution_index = src->getMaxBundleIndex();

    for (int i = min_distribution_index[0] ; i < max_distribution_index[0] ; ++ i) {
        for (int j = min_distribution_index[1] ; j < max_distribution_index[1] ; ++ j) {
            const int ci = (i - min_distribution_index[0]) * static_cast<int>(chunk_step);
            const int cj = (j - min_distribution_index[1]) * static_cast<int>(chunk_step);

            for (int k = 0 ; k < chunk_step ; ++k) {
                for (int l = 0 ; l < chunk_step ; ++l) {
                    const cslibs_math_2d::Point2d p(i * bundle_resolution + k * sampling_resolution,
                                                    j * bundle_resolution + l * sampling_resolution);

                    dst->at(static_cast<std::size_t>(ci + k),
                            static_cast<std::size_t>(cj + l)) =
                            src->sampleNonNormalized(
                                src->getInitialOrigin() * p, {{static_cast<int>(i), static_cast<int>(j)}});
                }
            }
        }
    }
}

inline void from(
        const cslibs_ndt_2d::static_maps::flat::Gridmap::Ptr  &src,
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr &dst,
        const double sampling_resolution)
{
    if (!src)
        return;

    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            src->getHeight() / sampling_resolution,
                            src->getWidth()  / sampling_resolution));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const std::size_t height = dst->getHeight();
    const std::size_t width  = dst->getWidth();
    const cslibs_ndt_2d::static_maps::flat::Gridmap::index_t min_index = src->getMinIndex();

    //// Todo: Check that code here!

    std::cerr << "converting " << width << " " << height << "\n";

    for(std::size_t i = 0 ; i < height ; ++i) {
      for(std::size_t j = 0 ; j < width ; ++j) {
        const cslibs_math_2d::Point2d  p(j * sampling_resolution,
                                         i * sampling_resolution);

        const cslibs_ndt_2d::static_maps::flat::Gridmap::index_t idx = {{min_index[0] + static_cast<int>(p(0) / src->getResolution()),
                                                                         min_index[1] + static_cast<int>(p(1) / src->getResolution())}};

        const double v = src->sampleNonNormalized(src->getOrigin() * p,idx);
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

    using dst_map_t = cslibs_gridmaps::static_maps::ProbabilityGridmap;
    dst.reset(new dst_map_t(src->getOrigin(),
                            sampling_resolution,
                            src->getHeight() / sampling_resolution,
                            src->getWidth()  / sampling_resolution));
    std::fill(dst->getData().begin(), dst->getData().end(), 0);

    const double bundle_resolution = src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinDistributionIndex();
    const index_t max_distribution_index = src->getMaxDistributionIndex();

    for (int i = min_distribution_index[0] ; i < max_distribution_index[0] ; ++ i) {
        for (int j = min_distribution_index[1] ; j < max_distribution_index[1] ; ++ j) {
            const int ci = (i - min_distribution_index[0]) * static_cast<int>(chunk_step);
            const int cj = (j - min_distribution_index[1]) * static_cast<int>(chunk_step);

            for (int k = 0 ; k < chunk_step ; ++k) {
                for (int l = 0 ; l < chunk_step ; ++l) {
                    const cslibs_math_2d::Point2d p(i * bundle_resolution + k * sampling_resolution,
                                                    j * bundle_resolution + l * sampling_resolution);

                    dst->at(static_cast<std::size_t>(ci + k),
                            static_cast<std::size_t>(cj + l)) =
                            src->sampleNonNormalized(
                                src->getInitialOrigin() * p, {{static_cast<int>(i), static_cast<int>(j)}}, inverse_model);
                }
            }
        }
    }
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_PROBABILITY_GRIDMAP_HPP
