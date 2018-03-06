#ifndef CSLIBS_NDT_3D_CONVERSION_POINTCLOUD_HPP
#define CSLIBS_NDT_3D_CONVERSION_POINTCLOUD_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace cslibs_ndt_3d {
namespace conversion {
inline void from(
        const cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr &src,
        pcl::PointCloud<pcl::PointXYZI>::Ptr &dst)
{
    if (!src)
        return;

    using dst_map_t = pcl::PointCloud<pcl::PointXYZI>;
    dst.reset(new dst_map_t());

    using index_t = std::array<int, 3>;
    const index_t min_distribution_index = src->getMinDistributionIndex();
    const index_t max_distribution_index = src->getMaxDistributionIndex();

    for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
        for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
            for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz) {
                const index_t bi({idx, idy, idz});
                if (const auto &b = src->getDistributionBundle(bi)) {
                    cslibs_math::statistics::Distribution<3, 3> d;
                    for (std::size_t i = 0 ; i < 8 ; ++i)
                        d += b->at(i)->getHandle()->data();
                    if (d.getN() == 0)
                        continue;

                    cslibs_math_3d::Point3d mean(d.getMean());
                    pcl::PointXYZI p;
                    p.x = static_cast<float>(mean(0));
                    p.y = static_cast<float>(mean(1));
                    p.z = static_cast<float>(mean(2));
                    p.intensity = static_cast<float>(src->sampleNonNormalized(mean));

                    dst->push_back(p);
                }
            }
        }
    }
}

inline void from(
        const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &src,
        pcl::PointCloud<pcl::PointXYZI>::Ptr &dst,
        const cslibs_gridmaps::utility::InverseModel::Ptr &ivm,
        const double &threshold = 0.169)
{
    if (!src)
        return;

    using dst_map_t = pcl::PointCloud<pcl::PointXYZI>;
    dst.reset(new dst_map_t());

    using index_t = std::array<int, 3>;
    const index_t min_distribution_index = src->getMinDistributionIndex();
    const index_t max_distribution_index = src->getMaxDistributionIndex();

    for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
        for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
            for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz) {
                const index_t bi({idx, idy, idz});
                if (const auto &b = src->getDistributionBundle(bi)) {
                    cslibs_math::statistics::Distribution<3, 3> d;
                    double occupancy = 0.0;

                    for (std::size_t i = 0 ; i < 8 ; ++i) {
                        occupancy += 0.125 * b->at(i)->getOccupancy(ivm);
                        if (b->at(i)->getDistribution())
                            d += *(b->at(i)->getDistribution());
                    }
                    if (d.getN() == 0 || occupancy < threshold)
                        continue;

                    cslibs_math_3d::Point3d mean(d.getMean());
                    pcl::PointXYZI p;
                    p.x = static_cast<float>(mean(0));
                    p.y = static_cast<float>(mean(1));
                    p.z = static_cast<float>(mean(2));
                    p.intensity = static_cast<float>(src->sampleNonNormalized(mean, ivm));

                    dst->push_back(p);
                }
            }
        }
    }
}
}
}

#endif // CSLIBS_NDT_3D_CONVERSION_POINTCLOUD_HPP
