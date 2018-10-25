#ifndef CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
#define CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP

#include <cslibs_ndt/matching/match.hpp>
#include <cslibs_ndt_3d/matching/gridmap_match_traits.hpp>

#include <cslibs_ndt_3d/matching/voxel.hpp>

namespace cslibs_ndt_3d {
namespace matching {
namespace dynamic_maps {
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                  const cslibs_ndt::matching::Parameter        &params,
                  double                                        resolution,
                  const cslibs_math_3d::Transform3d            &initial_transform,
                  cslibs_ndt::matching::Result<cslibs_math_3d::Transform3d> &r)
{
  using ndt_t = cslibs_ndt_3d::dynamic_maps::Gridmap;
  ndt_t ndt(ndt_t::pose_t(), resolution);
  ndt.insert(dst);
  r = cslibs_ndt::matching::match(src->begin(), src->end(), ndt, params, initial_transform);
}
}
}
}

#endif // CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
