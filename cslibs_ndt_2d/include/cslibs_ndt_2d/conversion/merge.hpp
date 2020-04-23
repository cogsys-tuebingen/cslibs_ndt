#ifndef CSLIBS_NDT_2D_FLATTEN_HPP
#define CSLIBS_NDT_2D_FLATTEN_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/mono_gridmap.hpp>

namespace cslibs_ndt_2d {
namespace conversion {
template <typename T>
inline typename cslibs_ndt_2d::static_maps::mono::Gridmap<T>::Ptr merge(
        const typename cslibs_ndt_2d::dynamic_maps::Gridmap<T>::Ptr& src)
{
    if (!src)
        return nullptr;

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinBundleIndex();
    const index_t max_distribution_index = src->getMaxBundleIndex();

    const std::array<std::size_t, 2> size = {{static_cast<std::size_t>(max_distribution_index[0] - min_distribution_index[0] + 1),
                                              static_cast<std::size_t>(max_distribution_index[1] - min_distribution_index[1] + 1)}};

    using src_map_t = cslibs_ndt_2d::dynamic_maps::Gridmap<T>;
    using dst_map_t = cslibs_ndt_2d::static_maps::mono::Gridmap<T>;

    typename dst_map_t::Ptr dst(new dst_map_t(src->getInitialOrigin(),
                                              src->getBundleResolution(),
                                              size,
                                              min_distribution_index));

    auto traverse = [&dst](const index_t &i, const typename src_map_t::distribution_bundle_t &b)
    {
        typename dst_map_t::distribution_t* dst_d = dst->getDistribution(i);
        for(const typename src_map_t::distribution_t *src_d : b) {
            *dst_d += *src_d;
        }
    };

    src->traverse(traverse);
    return dst;
}
}
}

#endif // CSLIBS_NDT_2D_FLATTEN_HPP
