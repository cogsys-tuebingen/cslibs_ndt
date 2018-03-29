#ifndef CSLIBS_NDT_2D_CONVERSION_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_CONVERSION_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/occupancy_gridmap.hpp>

namespace cslibs_ndt_2d {
namespace conversion {
inline cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr from(
        const cslibs_ndt_2d::static_maps::OccupancyGridmap::Ptr& src)
{
    if (!src)
        return nullptr;

    using src_map_t = cslibs_ndt_2d::static_maps::OccupancyGridmap;
    using dst_map_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap;
    typename dst_map_t::Ptr dst(new dst_map_t(src->getOrigin(),
                                              src->getResolution()));

    using index_t = std::array<int, 2>;
    src->traverse([&dst](const index_t &bi, const src_map_t::distribution_bundle_t &b){
        if (const typename dst_map_t::distribution_bundle_t* b_dst = dst->getDistributionBundle(bi)) {
            for (std::size_t i = 0 ; i < 4 ; ++ i)
                if (b.at(i)) {
                    const auto &handle = b.at(i)->getHandle();
                    if (handle->numFree() > 0 || handle->numOccupied() > 0)
                        *(b_dst->at(i)) = *handle;
                }
        }
    });

    return dst;
}

inline cslibs_ndt_2d::static_maps::OccupancyGridmap::Ptr from(
        const cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr& src)
{
    if (!src)
        return nullptr;

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinDistributionIndex();

    const std::array<std::size_t, 2> size =
    {{static_cast<std::size_t>(std::ceil((src->getMax()(0) - src->getMin()(0)) / src->getResolution())),
      static_cast<std::size_t>(std::ceil((src->getMax()(1) - src->getMin()(1)) / src->getResolution()))}};

    auto get_bundle_index = [&size, &min_distribution_index] (const index_t & bi) {
        return index_t{{bi[0] - min_distribution_index[0],
                        bi[1] - min_distribution_index[1]}};
    };

    using src_map_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap;
    using dst_map_t = cslibs_ndt_2d::static_maps::OccupancyGridmap;
    typename dst_map_t::Ptr dst(new dst_map_t(src->getOrigin(),
                                              src->getResolution(),
                                              size));

    src->traverse([&dst, &get_bundle_index](const index_t &bi, const src_map_t::distribution_bundle_t &b){
        const index_t bi_dst = get_bundle_index(bi);
        if (const typename dst_map_t::distribution_bundle_t* b_dst = dst->getDistributionBundle(bi_dst)) {
            for (std::size_t i = 0 ; i < 4 ; ++ i)
                if (b.at(i)) {
                    const auto &handle = b.at(i)->getHandle();
                    if (handle->numFree() > 0 || handle->numOccupied() > 0)
                        *(b_dst->at(i)) = *handle;
                }
        }
    });

    return dst;
}
}
}

#endif // CSLIBS_NDT_2D_CONVERSION_OCCUPANCY_GRIDMAP_HPP
