#ifndef CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP
#define CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_ndt_2d/conversion/gridmap.hpp>
#include <cslibs_ndt_2d/conversion/occupancy_gridmap.hpp>

#include <cslibs_gridmaps/static_maps/likelihood_field_gridmap.h>

namespace cslibs_ndt_2d {
namespace conversion {
inline cslibs_gridmaps::static_maps::LikelihoodFieldGridmap::Ptr from(
        const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr &src,
        const double &sampling_resolution)
{
    using src_map_t = cslibs_ndt_2d::static_maps::Gridmap;
    using dst_map_t = cslibs_gridmaps::static_maps::LikelihoodFieldGridmap;

    const typename src_map_t::Ptr static_src = cslibs_ndt_2d::conversion::from(src);
    typename dst_map_t::Ptr dst(new dst_map_t(static_src->getOrigin(),
                                              sampling_resolution,
                                              static_src->getHeight() / sampling_resolution,
                                              static_src->getWidth()  / sampling_resolution,
                                              0.0,   // ignored
                                              0.0)); // ignored

    const double bundle_resolution = static_src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);
    for (std::size_t i = 0 ; i < static_src->getBundleSize()[0] ; ++ i) {
        for (std::size_t j = 0 ; j < static_src->getBundleSize()[1] ; ++ j) {
            const int ci = i * chunk_step;
            const int cj = j * chunk_step;

            for (int k = 0 ; k < chunk_step ; ++k) {
                for (int l = 0 ; l < chunk_step ; ++l) {
                    const cslibs_math_2d::Point2d p(i * bundle_resolution + k * sampling_resolution,
                                                    j * bundle_resolution + l * sampling_resolution);

                    dst->at(static_cast<std::size_t>(ci + k),
                            static_cast<std::size_t>(cj + l)) =
                            static_src->sampleNonNormalized(
                                dst->getOrigin() * p, {{static_cast<int>(i), static_cast<int>(j)}});
                }
            }
        }
    }

    return dst;
}

inline cslibs_gridmaps::static_maps::LikelihoodFieldGridmap::Ptr from(
        const cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &src,
        const double &sampling_resolution,
        const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model)
{
    using src_map_t = cslibs_ndt_2d::static_maps::OccupancyGridmap;
    using dst_map_t = cslibs_gridmaps::static_maps::LikelihoodFieldGridmap;

    const typename src_map_t::Ptr static_src = cslibs_ndt_2d::conversion::from(src);
    typename dst_map_t::Ptr dst(new dst_map_t(static_src->getOrigin(),
                                              sampling_resolution,
                                              static_src->getHeight() / sampling_resolution,
                                              static_src->getWidth()  / sampling_resolution,
                                              0.0,   // ignored
                                              0.0)); // ignored

    const double bundle_resolution = static_src->getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);
    for (std::size_t i = 0 ; i < static_src->getBundleSize()[0] ; ++ i) {
        for (std::size_t j = 0 ; j < static_src->getBundleSize()[1] ; ++ j) {
            const int ci = i * chunk_step;
            const int cj = j * chunk_step;

            for (int k = 0 ; k < chunk_step ; ++k) {
                for (int l = 0 ; l < chunk_step ; ++l) {
                    const cslibs_math_2d::Point2d p(i * bundle_resolution + k * sampling_resolution,
                                                    j * bundle_resolution + l * sampling_resolution);

                    dst->at(static_cast<std::size_t>(ci + k),
                            static_cast<std::size_t>(cj + l)) =
                            static_src->sampleNonNormalized(
                                dst->getOrigin() * p, {{static_cast<int>(i), static_cast<int>(j)}}, inverse_model);
                }
            }
        }
    }

    return dst;
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_LIKELIHOOD_FIELD_GRIDMAP_HPP
