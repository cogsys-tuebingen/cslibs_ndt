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
        pcl::PointCloud<pcl::PointXYZI>::Ptr &dst,
        const bool &traversal = false)
{
    if (!src)
        return;

    using dst_map_t = pcl::PointCloud<pcl::PointXYZI>;
    dst.reset(new dst_map_t());

    using index_t = std::array<int, 3>;
    auto process_bundle = [&src, &dst](const index_t &bi) {
        if (const auto &b = src->getDistributionBundle(bi)) {
            cslibs_math::statistics::Distribution<3, 3> d;
            for (std::size_t i = 0 ; i < 8 ; ++i)
                d += b->at(i)->getHandle()->data();
            if (d.getN() == 0)
                return;

            cslibs_math_3d::Point3d mean(d.getMean());
            pcl::PointXYZI p;
            p.x = static_cast<float>(mean(0));
            p.y = static_cast<float>(mean(1));
            p.z = static_cast<float>(mean(2));
            p.intensity = static_cast<float>(src->sampleNonNormalized(mean));

            dst->push_back(p);
        }
    };

    if (traversal) {
        std::vector<index_t> indices;
        src->getBundleIndices(indices);
        for (auto &bi : indices)
            process_bundle(bi);
    } else {
        const index_t min_distribution_index = src->getMinDistributionIndex();
        const index_t max_distribution_index = src->getMaxDistributionIndex();
        if (min_distribution_index[0] == std::numeric_limits<int>::max() ||
                min_distribution_index[1] == std::numeric_limits<int>::max() ||
                min_distribution_index[2] == std::numeric_limits<int>::max() ||
                max_distribution_index[0] == std::numeric_limits<int>::min() ||
                min_distribution_index[1] == std::numeric_limits<int>::min() ||
                min_distribution_index[2] == std::numeric_limits<int>::min()) {
            dst = nullptr;
            return;
        }

        for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx)
            for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy)
                for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz)
                    process_bundle({{idx, idy, idz}});
    }
}

inline void from(
        const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &src,
        pcl::PointCloud<pcl::PointXYZI>::Ptr &dst,
        const cslibs_gridmaps::utility::InverseModel::Ptr &ivm,
        const bool &traversal = false,
        const double &threshold = 0.169)
{
    if (!src)
        return;

    using dst_map_t = pcl::PointCloud<pcl::PointXYZI>;
    dst.reset(new dst_map_t());

    using index_t = std::array<int, 3>;
    auto process_bundle = [&src, &dst, &ivm, &threshold](const index_t &bi) {
        if (const auto &b = src->getDistributionBundle(bi)) {
            cslibs_math::statistics::Distribution<3, 3> d;
            double occupancy = 0.0;

            for (std::size_t i = 0 ; i < 8 ; ++i) {
                occupancy += 0.125 * b->at(i)->getOccupancy(ivm);
                if (b->at(i)->getDistribution())
                    d += *(b->at(i)->getDistribution());
            }
            if (d.getN() == 0 || occupancy < threshold)
                return;

            cslibs_math_3d::Point3d mean(d.getMean());
            pcl::PointXYZI p;
            p.x = static_cast<float>(mean(0));
            p.y = static_cast<float>(mean(1));
            p.z = static_cast<float>(mean(2));
            p.intensity = static_cast<float>(src->sampleNonNormalized(mean, ivm));

            dst->push_back(p);
        }
    };

    if (traversal) {
        std::vector<index_t> indices;
        src->getBundleIndices(indices);
        for (auto &bi : indices)
            process_bundle(bi);
    } else {
        const index_t min_distribution_index = src->getMinDistributionIndex();
        const index_t max_distribution_index = src->getMaxDistributionIndex();
        if (min_distribution_index[0] == std::numeric_limits<int>::max() ||
                min_distribution_index[1] == std::numeric_limits<int>::max() ||
                min_distribution_index[2] == std::numeric_limits<int>::max() ||
                max_distribution_index[0] == std::numeric_limits<int>::min() ||
                min_distribution_index[1] == std::numeric_limits<int>::min() ||
                min_distribution_index[2] == std::numeric_limits<int>::min()) {
            dst = nullptr;
            return;
        }

        for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx)
            for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy)
                for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz)
                    process_bundle({{idx, idy, idz}});
    }
}
}
}

#endif // CSLIBS_NDT_3D_CONVERSION_POINTCLOUD_HPP
