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

template<typename T,
         typename ndt_t,
         typename = typename std::enable_if<std::is_same<ndt_t, cslibs_ndt_3d::dynamic_maps::Gridmap<T>>::value
                                            || std::is_same<ndt_t, cslibs_ndt_3d::static_maps::Gridmap<T>>::value>::type>
inline void rgbFrom(
        ndt_t &src,
        sensor_msgs::PointCloud2 &dst,
        const T sampling_resolution,
        const T& threshold = 0.196,
        const bool& allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using index_t = std::array<int, 3>;
    using point_t = typename ndt_t::point_t;
    using distribution_t = typename ndt_t::distribution_t;
    using distribution_bundle_t = typename ndt_t::distribution_bundle_t;
    auto sample = [](const distribution_t *d,
                     const point_t &p) -> T {
        return d ? d->data().sampleNonNormalized(p) : 0.0;
    };
    auto sample_bundle = [&sample](const distribution_bundle_t &b,
                                   const point_t &p) -> T {
        return 0.125 * (sample(b.at(0), p) +
                        sample(b.at(1), p) +
                        sample(b.at(2), p) +
                        sample(b.at(3), p) +
                        sample(b.at(4), p) +
                        sample(b.at(5), p) +
                        sample(b.at(6), p) +
                        sample(b.at(7), p));
    };

    typename cslibs_math_3d::PointcloudRGB3<T>::Ptr cloud(new typename cslibs_math_3d::PointcloudRGB3<T>);
    const T bundle_resolution = src.getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto process_bundle = [&src, &cloud, &sample_bundle, &chunk_step, &bundle_resolution, &sampling_resolution, &threshold](
            const index_t &bi, const distribution_bundle_t &b) {
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                for (int m = 0 ; m < chunk_step ; ++ m) {
                    const cslibs_math_3d::Point3<T> p(static_cast<T>(bi[0]) * bundle_resolution + (static_cast<T>(k)+0.5) * sampling_resolution,
                                                      static_cast<T>(bi[1]) * bundle_resolution + (static_cast<T>(l)+0.5) * sampling_resolution,
                                                      static_cast<T>(bi[2]) * bundle_resolution + (static_cast<T>(m)+0.5) * sampling_resolution);
                    if (sample_bundle(b, p) >= threshold)
                        cloud->insert(cslibs_math_3d::PointRGB3<T>(p));
                }
            }
        }
    };
    src.traverse(process_bundle);

    const T min_z = cloud->min().getPoint()(2);
    const T max_z = cloud->max().getPoint()(2);
    cloud->transform(ColorTransform<T>(min_z, max_z));
    cslibs_math_ros::sensor_msgs::conversion_3d::from<T>(cloud, dst);
}

template <typename T>
inline void rgbFrom(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const T sampling_resolution,
        const T& threshold = 0.196,
        const bool& allocate_all = false)
{
    if (!src)
        return;

    rgbFrom<T>(*src, dst, sampling_resolution, threshold, allocate_all);
}

template<typename T,
         typename ndt_t,
         typename = typename std::enable_if<std::is_same<ndt_t, cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>>::value
                                            || std::is_same<ndt_t, cslibs_ndt_3d::static_maps::OccupancyGridmap<T>>::value>::type>
inline void rgbFrom(
        ndt_t &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const T sampling_resolution,
        const T& threshold = 0.196,
        const bool& allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using index_t = std::array<int, 3>;
    using point_t = typename ndt_t::point_t;
    using distribution_t = typename ndt_t::distribution_t;
    using distribution_bundle_t = typename ndt_t::distribution_bundle_t;
    auto sample = [&ivm](const distribution_t *d,
                         const point_t &p) -> T {
        auto evaluate = [&ivm, d, p] {
            const auto &handle = d;
            return d && handle->getDistribution() ?
                        handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : 0.0;
        };
        return d ? evaluate() : 0.0;
    };
    auto sample_bundle = [&sample](const distribution_bundle_t &b,
                                   const point_t &p) -> T {
        return 0.125 * (sample(b.at(0), p) +
                        sample(b.at(1), p) +
                        sample(b.at(2), p) +
                        sample(b.at(3), p) +
                        sample(b.at(4), p) +
                        sample(b.at(5), p) +
                        sample(b.at(6), p) +
                        sample(b.at(7), p));
    };


    typename cslibs_math_3d::PointcloudRGB3<T>::Ptr cloud(new typename cslibs_math_3d::PointcloudRGB3<T>);
    const T bundle_resolution = src.getBundleResolution();
    const int chunk_step = static_cast<int>(bundle_resolution / sampling_resolution);

    auto process_bundle = [&src, &cloud, &sample_bundle, &chunk_step, &bundle_resolution, &sampling_resolution, &threshold](
            const index_t &bi, const distribution_bundle_t &b) {
        for (int k = 0 ; k < chunk_step ; ++ k) {
            for (int l = 0 ; l < chunk_step ; ++ l) {
                for (int m = 0 ; m < chunk_step ; ++ m) {
                    const cslibs_math_3d::Point3<T> p(static_cast<T>(bi[0]) * bundle_resolution + (static_cast<T>(k)+0.5) * sampling_resolution,
                                                      static_cast<T>(bi[1]) * bundle_resolution + (static_cast<T>(l)+0.5) * sampling_resolution,
                                                      static_cast<T>(bi[2]) * bundle_resolution + (static_cast<T>(m)+0.5) * sampling_resolution);
                    if (sample_bundle(b, p) >= threshold)
                        cloud->insert(cslibs_math_3d::PointRGB3<T>(p));
                }
            }
        }
    };
    src.traverse(process_bundle);

    const T min_z = cloud->min().getPoint()(2);
    const T max_z = cloud->max().getPoint()(2);
    cloud->transform(ColorTransform<T>(min_z, max_z));
    cslibs_math_ros::sensor_msgs::conversion_3d::from<T>(cloud, dst);
}

template <typename T>
inline void rgbFrom(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const T sampling_resolution,
        const T& threshold = 0.196,
        const bool& allocate_all = false)
{
    if (!src)
        return;

    rgbFrom<T>(*src, dst, ivm, sampling_resolution, threshold, allocate_all);
}

}
}

#endif // CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_RGB_HPP
