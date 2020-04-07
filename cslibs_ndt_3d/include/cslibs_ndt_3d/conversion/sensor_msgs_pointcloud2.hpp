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

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void from(
        const cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t> &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const bool &allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using point_t = typename ndt_t::point_t;
    using distribution_t = typename ndt_t::distribution_t;
    auto sample = [](const distribution_t *d,
                     const point_t &p) -> T {
        return d ? d->sampleNonNormalized(p) : 0.0;
        //return d && d->getDistribution() ? d->getDistribution()->sampleNonNormalized(p) : 0.0;
    };

    std::vector<float> tmp;
    const auto& origin = transform * src.getInitialOrigin();
    auto process_item = [&src, &tmp, &origin, &sample](const index_t &bi, const distribution_t& d) {
        /*const auto dd = d.getDistribution();
        if (!dd)
            return;*/

        cslibs_math_3d::Point3<T> mean(/*dd->*/d.getMean());
        cslibs_math_3d::Point3<T> p = origin * mean;
        tmp.emplace_back(static_cast<float>(p(0)));
        tmp.emplace_back(static_cast<float>(p(1)));
        tmp.emplace_back(static_cast<float>(p(2)));
        tmp.emplace_back(static_cast<float>(sample(&d,mean)));
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }
    from(tmp, dst);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const bool &allocate_all = false)
{
    if (!src)
        return;

    from(*src, dst, transform, allocate_all);
}

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void from(
        const cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t> &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T &threshold = 0.169,
        const bool &allocate_all = false)
{
    if (allocate_all)
        src.allocatePartiallyAllocatedBundles();

    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using point_t = typename ndt_t::point_t;
    using distribution_t = typename ndt_t::distribution_t;
    auto sample = [&ivm](const distribution_t *d,
                         const point_t &p) -> T {
        auto evaluate = [&ivm, &d, &p] {
            const auto &handle = d;
            return handle->getDistribution() ?
                        handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : T(0.0);
        };
        return d ? evaluate() : T(0.0);
    };

    std::vector<float> tmp;
    const auto& origin = transform * src.getInitialOrigin();
    auto process_item = [&src, &tmp, &origin, &sample, &ivm, &threshold](const index_t &bi, const distribution_t& d) {
        if (d.getOccupancy(ivm) < threshold)
            return;

        const auto dd = d.getDistribution();
        if (!dd)
            return;
        cslibs_math_3d::Point3<T> mean(dd->getMean());
        cslibs_math_3d::Point3<T> p = origin * mean;
        tmp.emplace_back(static_cast<float>(p(0)));
        tmp.emplace_back(static_cast<float>(p(1)));
        tmp.emplace_back(static_cast<float>(p(2)));
        tmp.emplace_back(static_cast<float>(sample(&d,mean)));
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }
    from(tmp, dst);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        sensor_msgs::PointCloud2 &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T &threshold = 0.169,
        const bool &allocate_all = false)
{
    if (!src)
        return;

    from(*src, dst, ivm, transform, threshold, allocate_all);
}

}
}

#endif // CSLIBS_NDT_3D_CONVERSION_SENSOR_MSGS_POINTCLOUD2_HPP
