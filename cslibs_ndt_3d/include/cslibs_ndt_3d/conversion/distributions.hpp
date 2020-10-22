#ifndef CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
#define CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_math/color/color.hpp>
#include <cslibs_math/common/angle.hpp>
#include <visualization_msgs/MarkerArray.h>

namespace cslibs_ndt_3d {
namespace conversion {

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void from(
        const cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t> &src,
        visualization_msgs::MarkerArray &dst,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
{
    using src_map_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using point_t = typename src_map_t::point_t;
    using distribution_t = typename src_map_t::distribution_t;

    visualization_msgs::Marker marker;
    marker.header.stamp = time;
    marker.header.frame_id = frame;
    marker.ns = "distributions";
    marker.action = visualization_msgs::Marker::DELETEALL;
    dst.markers.push_back(marker);

    marker.action = visualization_msgs::Marker::ADD;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.color.a = 1;
    marker.lifetime = ros::Duration(2000.);
    marker.id++;

    const auto& origin = transform * src.getInitialOrigin();
    const T min_height = (origin * src.getMin())(2);
    const T max_height = (origin * src.getMax())(2);

    auto process_item = [&dst,&marker,&origin,&min_height,&max_height](
            const index_t &bi, const distribution_t& d) {        
        /*const auto& data = d.getDistribution();
        if (!data) return;*/

        const auto& mean = /*data->*/d.getMean();
        const auto& p = origin * point_t(mean);

        marker.pose.position.x = p(0);
        marker.pose.position.y = p(1);
        marker.pose.position.z = p(2);

        typename distribution_t::distribution_t::eigen_values_t eval;
        typename distribution_t::distribution_t::eigen_vectors_t evec;
        if (!/*data->*/d.getEigenValuesVectors(eval,evec,true))
            return;

        //const auto& evec = data.getEigenVectors();
        const Eigen::Quaternion<T> orientation =
                origin.rotation().toEigen() * Eigen::Quaternion<T>(evec);

        marker.pose.orientation.x = orientation.x();
        marker.pose.orientation.y = orientation.y();
        marker.pose.orientation.z = orientation.z();
        marker.pose.orientation.w = orientation.w();

        //const auto& eval = data.getEigenValues();
        marker.scale.x = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(0)));
        marker.scale.y = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(1)));
        marker.scale.z = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(2)));

        cslibs_math::color::Color<T> color =
                cslibs_math::color::interpolateColor(p(2), min_height, max_height);
        marker.color.a = 1;
        marker.color.r = color.r;
        marker.color.g = color.g;
        marker.color.b = color.b;

        dst.markers.push_back(marker);
        marker.id++;
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        visualization_msgs::MarkerArray::Ptr &dst,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>())
{
    if (!src)
        return;
    dst.reset(new visualization_msgs::MarkerArray());

    from(*src,*dst,time,frame,transform);
}

template <cslibs_ndt::map::tags::option option_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
inline void from(
        const cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t> &src,
        visualization_msgs::MarkerArray &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T occupancy_threshold = 0.5)
{
    using src_map_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::OccupancyDistribution,T,backend_t>;
    using index_t = std::array<int, 3>;
    using point_t = typename src_map_t::point_t;
    using distribution_t = typename src_map_t::distribution_t;

    visualization_msgs::Marker marker;
    marker.header.stamp = time;
    marker.header.frame_id = frame;
    marker.ns = "distributions";
    marker.action = visualization_msgs::Marker::DELETEALL;
    dst.markers.push_back(marker);

    marker.action = visualization_msgs::Marker::ADD;
    marker.type = visualization_msgs::Marker::SPHERE;
    marker.color.a = 1;
    marker.lifetime = ros::Duration(2000.);
    marker.id++;

    const auto& origin = transform * src.getInitialOrigin();
    const T min_height = (origin * src.getMin())(2);
    const T max_height = (origin * src.getMax())(2);

    auto process_item = [&dst,&ivm,&occupancy_threshold,&marker,&origin,&min_height,&max_height](
            const index_t &bi, const distribution_t& d) {
        const T occupancy = d.getOccupancy(ivm);
        if (occupancy < occupancy_threshold)
            return;

        const auto& data = d.getDistribution();
        if (!data) return;

        const auto& mean = data->getMean();
        const auto& p = origin * point_t(mean);

        marker.pose.position.x = p(0);
        marker.pose.position.y = p(1);
        marker.pose.position.z = p(2);

        typename distribution_t::distribution_t::eigen_values_t eval;
        typename distribution_t::distribution_t::eigen_vectors_t evec;
        if (!data->getEigenValuesVectors(eval,evec,true))
            return;

        //const auto& evec = data->getEigenVectors();
        const Eigen::Quaternion<T> orientation =
                origin.rotation().toEigen() * Eigen::Quaternion<T>(evec);

        marker.pose.orientation.x = orientation.x();
        marker.pose.orientation.y = orientation.y();
        marker.pose.orientation.z = orientation.z();
        marker.pose.orientation.w = orientation.w();

        //const auto& eval = data->getEigenValues();
        marker.scale.x = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(0)));
        marker.scale.y = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(1)));
        marker.scale.z = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(2)));

        cslibs_math::color::Color<T> color =
                cslibs_math::color::interpolateColor(p(2), min_height, max_height);
        marker.color.a = occupancy;
        marker.color.r = color.r;
        marker.color.g = color.g;
        marker.color.b = color.b;

        dst.markers.push_back(marker);
        marker.id++;
    };

    const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        visualization_msgs::MarkerArray::Ptr &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr& ivm,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T occupancy_threshold = 0.5)
{
    if (!src || !ivm)
        return;
    dst.reset(new visualization_msgs::MarkerArray());

    from(*src,*dst,ivm,time,frame,transform,occupancy_threshold);
}

}
}

#endif // CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
