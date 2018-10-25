#ifndef CSLIBS_NDT_2D_FLATTEN_HPP
#define CSLIBS_NDT_2D_FLATTEN_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/mono_gridmap.hpp>

namespace cslibs_ndt_2d {
namespace conversion {
inline cslibs_ndt_2d::static_maps::mono::Gridmap::Ptr merge(
        const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr& src)
{
    if (!src)
        return nullptr;

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinBundleIndex();
    const index_t max_distribution_index = src->getMaxBundleIndex();

    const std::array<std::size_t, 2> size = {{static_cast<std::size_t>(max_distribution_index[0] - min_distribution_index[0] + 1),
                                              static_cast<std::size_t>(max_distribution_index[1] - min_distribution_index[1] + 1)}};

    using src_map_t = cslibs_ndt_2d::dynamic_maps::Gridmap;
    using dst_map_t = cslibs_ndt_2d::static_maps::mono::Gridmap;

    typename dst_map_t::Ptr dst(new dst_map_t(src->getInitialOrigin(),
                                              src->getBundleResolution(),
                                              size,
                                              min_distribution_index));

    auto traverse = [&min_distribution_index, &dst](const index_t &i, const src_map_t::distribution_bundle_t &b)
    {
        dst_map_t::distribution_t* dst_d = dst->getDistribution(i);
        for(const src_map_t::distribution_t *src_d : b) {
            dst_d->data() += src_d->data();
        }
    };

    src->traverse(traverse);
    return dst;
}
}
}


#endif // CSLIBS_NDT_2D_FLATTEN_HPP
