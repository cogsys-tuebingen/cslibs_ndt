#ifndef CSLIBS_NDT_2D_FLATTEN_HPP
#define CSLIBS_NDT_2D_FLATTEN_HPP

#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/static_maps/flat_gridmap.hpp>

namespace cslibs_ndt_2d {
namespace conversion {
inline cslibs_ndt_2d::static_maps::flat::Gridmap::Ptr flatten(
        const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr& src)
{
    if (!src)
        return nullptr;

    using index_t = std::array<int, 2>;
    const index_t min_distribution_index = src->getMinDistributionIndex();

    const cslibs_math_2d::Point2d min_corner = src->getMin();
    const cslibs_math_2d::Point2d max_corner = src->getMax();
    const double resolution = src->getBundleResolution();


    const std::array<std::size_t, 2> size =
    {{static_cast<std::size_t>(std::ceil((max_corner(0) - min_corner(0)) / resolution)),
      static_cast<std::size_t>(std::ceil((max_corner(1) - min_corner(1)) / resolution))}};

    using src_map_t = cslibs_ndt_2d::dynamic_maps::Gridmap;
    using dst_map_t = cslibs_ndt_2d::static_maps::flat::Gridmap;

    typename dst_map_t::Ptr dst(new dst_map_t(src->getOrigin(),
                                              src->getResolution(),
                                              size));
    auto traverse = [&min_distribution_index, &dst](const index_t &i, const src_map_t::distribution_bundle_t &b)
    {
      const index_t index = {{i[0] - min_distribution_index[0], i[1] - min_distribution_index[1]}};
      dst_map_t::distribution_t* dst_d = dst->getDistribution(index);
      for(const src_map_t::distribution_t *src_d : b) {
        dst_d->getHandle()->data() += src_d->getHandle()->data();
      }

    };
    src->traverse(traverse);
    return dst;
}
}
}


#endif // CSLIBS_NDT_2D_FLATTEN_HPP
