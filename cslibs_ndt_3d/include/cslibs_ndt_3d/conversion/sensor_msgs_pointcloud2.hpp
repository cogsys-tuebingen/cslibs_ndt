#ifndef CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_HPP
#define CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/gridmap.hpp>
#include <cslibs_ndt_3d/static_maps/occupancy_gridmap.hpp>

#include <sensor_msgs/PointCloud2.h>

namespace cslibs_ndt_3d {
namespace conversion {
inline void from(
        const std::vector<float> &tmp,
        sensor_msgs::PointCloud2 &dst)
{
    // metadata
    dst.width        = tmp.size() / 4;
    dst.height       = 1;
    dst.is_dense     = false;
    dst.is_bigendian = false;
    dst.point_step   = 4 * sizeof(float);
    dst.row_step     = static_cast<uint32_t>(sizeof(float) * tmp.size());

    // fields x y z intensity
    dst.fields.resize(4);
    dst.fields[0].name = "x";
    dst.fields[1].name = "y";
    dst.fields[2].name = "z";
    dst.fields[3].name = "intensity";
    for (int i = 0; i < 4; ++i) {
        dst.fields[i].offset   = i * sizeof(float);
        dst.fields[i].datatype = sensor_msgs::PointField::FLOAT32;
        dst.fields[i].count    = dst.width;
    }

    // data
    std::size_t data_size = sizeof(float) * tmp.size();
    dst.data.resize(data_size);
    memcpy(&dst.data[0], &tmp[0], data_size);
}

template<typename T,
         typename ndt_t,
         typename = typename std::enable_if<std::is_same<ndt_t, cslibs_ndt_3d::dynamic_maps::Gridmap<T>>::value
                                            || std::is_same<ndt_t, cslibs_ndt_3d::static_maps::Gridmap<T>>::value>::type>
inline void from(
        ndt_t &src,
        sensor_msgs::PointCloud2 &dst,
        const bool &allocate_all = true,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
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

    std::vector<float> tmp;
    auto process_bundle = [&src, &tmp, &transform, &sample_bundle](const index_t &bi, const distribution_bundle_t &b) {
        typename distribution_t::distribution_t d;
        for (std::size_t i = 0 ; i < 8 ; ++i)
            d += b.at(i)->data();
        if (d.getN() == 0)
            return;

        cslibs_math_3d::Point3<T> mean(d.getMean());
        cslibs_math_3d::Point3<T> p = transform * mean;
        tmp.emplace_back(static_cast<float>(p(0)));
        tmp.emplace_back(static_cast<float>(p(1)));
        tmp.emplace_back(static_cast<float>(p(2)));
        tmp.emplace_back(static_cast<float>(sample_bundle(b, mean)));
    };
    src.traverse(process_bundle);
    from(tmp, dst);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const bool &allocate_all = true,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
{
    if (!src)
        return;

    from<T>(*src, dst, allocate_all, transform);
}

template<typename T,
         typename ndt_t,
         typename = typename std::enable_if<std::is_same<ndt_t, cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>>::value
                                            || std::is_same<ndt_t, cslibs_ndt_3d::static_maps::OccupancyGridmap<T>>::value>::type>
inline void from(
        ndt_t &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const T &threshold = 0.169,
        const bool &allocate_all = true,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
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

    std::vector<float> tmp;
    auto process_bundle = [&src, &tmp, &ivm, &threshold, &transform, &sample_bundle](const index_t &bi, const distribution_bundle_t &b) {
        typename distribution_t::distribution_t d;//cslibs_math::statistics::Distribution<T, 3, 3> d;
        T occupancy = 0.0;

        for (std::size_t i = 0 ; i < 8 ; ++i) {
            const auto &handle = b.at(i);
            occupancy += 0.125 * handle->getOccupancy(ivm);
            if (const auto &d_tmp = handle->getDistribution())
                d += *d_tmp;
        }
        if (d.getN() == 0 || occupancy < threshold)
            return;

        cslibs_math_3d::Point3<T> mean(d.getMean());
        cslibs_math_3d::Point3<T> p = transform * mean;
        tmp.emplace_back(static_cast<float>(p(0)));
        tmp.emplace_back(static_cast<float>(p(1)));
        tmp.emplace_back(static_cast<float>(p(2)));
        tmp.emplace_back(static_cast<float>(sample_bundle(b, mean)));
    };
    src.traverse(process_bundle);
    from(tmp, dst);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const T &threshold = 0.169,
        const bool &allocate_all = true,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
{
    if (!src)
        return;

    from<T>(*src, dst, ivm, threshold, allocate_all, transform);
}

}
}

#endif // CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_HPP
