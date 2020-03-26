#ifndef CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_RGB_HPP
#define CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_RGB_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/occupancy_gridmap.hpp>
#include <cslibs_math_ros/sensor_msgs/conversion_3d.hpp>

#include <sensor_msgs/PointCloud2.h>

namespace cslibs_ndt_3d {
namespace conversion {
template <typename T>
struct ColorTransform {
    ColorTransform(const double& min_z, const double& max_z) : min_z_(min_z), max_z_(max_z) {}
    cslibs_math_3d::PointRGB3<T> apply(const cslibs_math_3d::PointRGB3<T>& p) const
    {
        return cslibs_math_3d::PointRGB3<T>(p.getPoint(), 1.0f, cslibs_math::color::interpolateColor<T>(p.getPoint()(2), min_z_, max_z_));
    }

private:
    const double min_z_;
    const double max_z_;
};

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void rgbFrom(
        cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t> &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const bool& allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using distribution_t = typename ndt_t::distribution_t;

    cslibs_math_3d::PointcloudRGB3<T> cloud;
    const auto& origin = transform * src.getInitialOrigin();
    auto process_item = [&src, &cloud, &origin](const index_t &bi, const distribution_t& d) {
        cslibs_math_3d::Point3<T> mean(d.data().getMean());
        cloud.insert(cslibs_math_3d::PointRGB3<T>(origin * mean));
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }

    const T min_z = cloud.min().getPoint()(2);
    const T max_z = cloud.max().getPoint()(2);
    cloud.transform(ColorTransform<T>(min_z, max_z));
    cslibs_math_ros::sensor_msgs::conversion_3d::from<T>(cloud, dst);
}

template <typename T>
inline void rgbFrom(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const bool& allocate_all = false)
{
    if (!src)
        return;

    rgbFrom(*src, dst, transform, allocate_all);
}

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void rgbFrom(
        cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t> &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T& threshold = 0.169,
        const bool& allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using distribution_t = typename ndt_t::distribution_t;

    cslibs_math_3d::PointcloudRGB3<T> cloud;
    const auto& origin = transform * src.getInitialOrigin();
    auto process_item = [&src, &ivm, &cloud, &origin, &threshold](const index_t &bi, const distribution_t& d) {
        const auto dd = d.getDistribution();
        if (!dd || d.getOccupancy(ivm) < threshold)
            return;
        cslibs_math_3d::PointRGB3<T> mean(dd->getMean());
        cloud.insert(cslibs_math_3d::PointRGB3<T>(origin * mean));
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }

    const T min_z = cloud.min().getPoint()(2);
    const T max_z = cloud.max().getPoint()(2);
    cloud.transform(ColorTransform<T>(min_z, max_z));
    cslibs_math_ros::sensor_msgs::conversion_3d::from<T>(cloud, dst);
}

template <typename T>
inline void rgbFrom(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T& threshold = 0.169,
        const bool& allocate_all = false)
{
    if (!src)
        return;

    rgbFrom<T>(*src, dst, ivm, transform, threshold, allocate_all);
}

}
}

#endif // CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_RGB_HPP
