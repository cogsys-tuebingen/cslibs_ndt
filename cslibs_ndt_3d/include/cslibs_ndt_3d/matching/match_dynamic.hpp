#ifndef CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
#define CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP

#include <cslibs_ndt/matching/match.hpp>
#include <cslibs_ndt/matching/voxel.hpp>

#include <cslibs_ndt_3d/matching/gridmap_match_traits.hpp>
#include <cslibs_ndt_3d/matching/icp_params.hpp>
#include <cslibs_ndt_3d/matching/icp_result.hpp>
#include <cslibs_ndt_3d/matching/icp.hpp>

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

inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr          &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr          &dst,
                  const cslibs_ndt_3d::matching::ParametersWithICP      &params,
                  double                                                 resolution,
                  const cslibs_math_3d::Transform3d                     &initial_transform,
                  cslibs_ndt_3d::matching::ResultWithICP                &r)
{
    using ndt_t         = cslibs_ndt_3d::dynamic_maps::Gridmap;
    using voxel_grid_t  = cslibs_ndt::matching::VoxelGrid<3>;
    using voxel_t       = cslibs_ndt::matching::Voxel<3>;

    auto create_voxeled_cloud = [resolution](const cslibs_math_3d::Pointcloud3d::ConstPtr &src)
    {
        const voxel_t::index_t min_index = voxel_t::getIndex(src->min(), resolution);
        const voxel_t::index_t max_index = voxel_t::getIndex(src->max(), resolution);
        const voxel_t::size_t  size      = {{static_cast<size_t>(max_index[0] - min_index[0] + 1),
                                             static_cast<size_t>(max_index[1] - min_index[1] + 1),
                                             static_cast<size_t>(max_index[2] - min_index[2] + 1)}};

        voxel_grid_t::Ptr voxel_grid(new voxel_grid_t::type);
        voxel_grid->set<cis::option::tags::array_offset>(min_index[0], min_index[1], min_index[2]);
        voxel_grid->set<cis::option::tags::array_size>(size[0], size[1], size[2]);

        const cslibs_math_3d::Pointcloud3d::points_t &pts = src->getPoints();
        for(const auto &p : pts) {
            voxel_grid->insert(voxel_t::getIndex(p, resolution), voxel_t(p));
        }

        cslibs_math_3d::Pointcloud3d::Ptr voxeled_cloud(new cslibs_math_3d::Pointcloud3d);
        auto traverse = [&voxeled_cloud](const voxel_t::index_t, voxel_t &voxel )
        {
                voxeled_cloud->insert(voxel.mean());
        };
        voxel_grid->traverse(traverse);
        return voxeled_cloud;

    };

    /// here we voxel the input clouds, to apply icp up front
    cslibs_ndt_3d::matching::impl::icp::apply(create_voxeled_cloud(src),
                                              create_voxeled_cloud(dst),
                                              params,
                                              initial_transform,
                                              r);


    ndt_t ndt(ndt_t::pose_t(), resolution);
    ndt.insert(dst);
    r.assign(cslibs_ndt::matching::match(src->begin(), src->end(), ndt, params, initial_transform));
}
}
}
}

#endif // CSLIBS_NDT_3D_MATCH_DYNAMIC_HPP
