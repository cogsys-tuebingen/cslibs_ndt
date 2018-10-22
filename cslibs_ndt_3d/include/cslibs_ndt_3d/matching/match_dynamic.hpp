#ifndef CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
#define CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP

#include <cslibs_ndt_3d/matching/match.hpp>
#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>

namespace cslibs_ndt_3d {
namespace matching {
namespace dynamic_maps {
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                  const Parameters                             &params,
                  Result                                       &r)
{
  using ndt_t = cslibs_ndt_3d::dynamic_maps::Gridmap;
  ndt_t::Ptr ndt(new ndt_t(ndt_t::pose_t(), params.resolution()));
  ndt->insert(dst);
  impl::match<ndt_t>(src, ndt, params, r);
}

template<std::size_t Ts>
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                  const Parameters                             &params,
                  Result                                       &r)
{
  using ndt_t = cslibs_ndt_3d::dynamic_maps::Gridmap;
  ndt_t::Ptr ndt(new ndt_t(ndt_t::pose_t(), params.resolution()));
  ndt->insert(dst);
  impl::match<ndt_t, Ts>(src, ndt, params, r);
}
}
}
}

#endif // CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
