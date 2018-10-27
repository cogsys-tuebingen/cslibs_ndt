#ifndef MATCH_STATIC_HPP
#define MATCH_STATIC_HPP

#include <cslibs_ndt/matching/match.hpp>
#include <cslibs_ndt_3d/matching/gridmap_match_traits.hpp>

#include <cslibs_ndt/matching/voxel.hpp>

namespace cslibs_ndt_3d {
namespace matching {
namespace static_maps {

inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                  const cslibs_ndt::matching::Parameter        &params,
                  double                                        resolution,
                  const cslibs_math_3d::Transform3d            &initial_transform,
                  cslibs_ndt::matching::Result<cslibs_math_3d::Transform3d> &r)
{
    using ndt_t   = ::cslibs_ndt_3d::static_maps::Gridmap;
    using size_t  = ndt_t::size_t;
    using index_t = ndt_t::index_t;

    const auto min = dst->min();
    const auto max = dst->max();
    const index_t min_index = {{static_cast<int>(min(0) / resolution),
                                static_cast<int>(min(1) / resolution),
                                static_cast<int>(min(2) / resolution)}};
    const size_t size       = {{static_cast<std::size_t>(static_cast<int>(max(0) / resolution) - min_index[0]) + 1,
                                static_cast<std::size_t>(static_cast<int>(max(1) / resolution) - min_index[1]) + 1,
                                static_cast<std::size_t>(static_cast<int>(max(2) / resolution) - min_index[2]) + 1}};

    ndt_t ndt(ndt_t::pose_t(), resolution, size, min_index);
    ndt.insert(dst);
    r = cslibs_ndt::matching::match(src->begin(), src->end(), ndt, params, initial_transform);
}
}
}
}

#endif // MATCH_STATIC_HPP
