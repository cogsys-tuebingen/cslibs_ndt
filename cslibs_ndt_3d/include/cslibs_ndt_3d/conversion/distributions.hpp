#ifndef CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
#define CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_math/color/color.hpp>
#include <cslibs_math/common/angle.hpp>
#include <visualization_msgs/MarkerArray.h>

namespace cslibs_ndt_3d {
namespace conversion {
/*template <template <typename, std::size_t, std::size_t> typename distribution_t, typename T>
inline Distribution from(const distribution_t<T, 3, 3> &d,
                         const int &id,
                         const T &prob)
{
    Distribution distr;
    distr.id.data = id;
    for (int i = 0; i < 3; ++ i) {
        distr.mean[i].data          = d.getMean()(i);
        distr.eigen_values[i].data  = d.getEigenValues()(i);
    }
    for (int i = 0; i < 9; ++ i) {
        distr.eigen_vectors[i].data = d.getEigenVectors()(i);
        distr.covariance[i].data    = d.getCovariance()(i);
    }
    distr.prob.data = prob;
    return distr;
}*/

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

template <typename T>
inline void from(
        const cslibs_ndt_3d::dynamic_maps::Gridmap<T> &src,
        visualization_msgs::MarkerArray &dst,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_math_3d::Pose3<T> &transform = typename cslibs_math_3d::Pose3<T>(),
        const T& alpha_step = 45.,
        const T& beta_step = 45.)
{
    const T alpha_step_rad   = cslibs_math::common::angle::toRad(alpha_step);
    const T beta_step_rad    = cslibs_math::common::angle::toRad(beta_step);
    static const T angle_max = cslibs_math::common::angle::toRad(360.);

    using src_map_t = cslibs_ndt_3d::dynamic_maps::Gridmap<T>;
    using index_t = std::array<int, 3>;
    using point_t = typename src_map_t::point_t;
    using pose_t = typename src_map_t::pose_t;
    using distribution_t = typename src_map_t::distribution_t;
    using bundle_t = typename src_map_t::distribution_bundle_t;

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

    auto process_item = [&marker,&dst,&alpha_step_rad,&beta_step_rad,&origin,&min_height,&max_height](
            const index_t &bi, const distribution_t& d) {
        const auto& data = d.data();
        const auto& mean = data.getMean();
        const auto& p = origin * point_t(mean);

        marker.pose.position.x = p(0);
        marker.pose.position.y = p(1);
        marker.pose.position.z = p(2);

        const auto& evec = data.getEigenVectors();
        const Eigen::Quaternion<T> orientation = origin.rotation().toEigen() * Eigen::Quaternion<T>(evec);

        marker.pose.orientation.x = orientation.x();
        marker.pose.orientation.y = orientation.y();
        marker.pose.orientation.z = orientation.z();
        marker.pose.orientation.w = orientation.w();

        const auto& eval = data.getEigenValues();
        marker.scale.x = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(0)));
        marker.scale.y = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(1)));
        marker.scale.z = std::max(static_cast<T>(1e-4),T(2.)*std::sqrt(eval(2)));

        cslibs_math::color::Color<T> color = cslibs_math::color::interpolateColor(p(2), min_height, max_height);
        marker.color.a = 1;
        marker.color.r = color.r;
        marker.color.g = color.g;
        marker.color.b = color.b;

        dst.markers.push_back(marker);
        marker.id++;
    };
    std::array<std::set<index_t>,src_map_t::bin_count> closed;
    auto process = [&process_item,&closed](const index_t& bi, const bundle_t b) {
        /*cslibs_ndt::utility::apply_indices<src_map_t::bin_count,3>(
                    bi, [&b,&closed,&process_item](const std::size_t& i, const index_t& index) {
            if (closed[i].find(index) == closed[i].end()) {
                if (b.at(i))
                    process_item(index,*(b.at(i)));
                closed[i].insert(index);
            }
        });/*/
        distribution_t d;
        for (std::size_t i=0; i<src_map_t::bin_count; ++i)
            if (b.at(i))
                d.data() += b.at(i)->data();
        process_item(bi,d);//*/
    };

    /*const auto& storages = src.getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }/*/
    src.traverse(process);//*/
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        visualization_msgs::MarkerArray::Ptr &dst,
        const ros::Time& time,
        const std::string &frame,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr& ivm)
{
    if (!src || !ivm)
        return;

    using src_map_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>;
    dst.reset(new visualization_msgs::MarkerArray());

    visualization_msgs::Marker marker;
    marker.header.stamp = time;
    marker.header.frame_id = frame;
    marker.action = visualization_msgs::Marker::DELETEALL;
    marker.ns = "distributions";
    dst->markers.push_back(marker);

    using index_t = std::array<int, 3>;
    marker.action = visualization_msgs::Marker::ADD;
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.scale.x = 0.02;
    marker.scale.y = 1;
    marker.scale.z = 1;
    marker.color.a = 1;
    marker.lifetime = ros::Duration();

    const auto& origin = src->getOrigin();
    const auto& min_height = (origin * src->getMin())(2);
    const auto& max_height = (origin * src->getMax())(2);

    auto to_point = [&origin](const Eigen::Matrix<T,3,1>& p) {
        geometry_msgs::Point p_res;
        const auto& pt = origin * p;
        p_res.x = pt(0);
        p_res.y = pt(1);
        p_res.z = pt(2);
        return p_res;
    };
    using distribution_t = typename src_map_t::distribution_t;
    auto process_item = [&dst,&marker,&to_point,&origin,&min_height,&max_height,&ivm](
            const index_t &bi, const distribution_t& d) {
        const auto& data = d.getDistribution();
        if (!data) return;

        const auto& mean = data->getMean();
        const auto& evec = data->getEigenVectors();
        const auto& eval = data->getEigenValues();

        const auto height = (origin * mean)(2);
        cslibs_math::color::Color<T> color = cslibs_math::color::interpolateColor(height, min_height, max_height);
        marker.color.a = d.getOccupancy(ivm);

        for (std::size_t i=0; i<3; ++i) {
            T scale = eval(i);
            if (scale < 1e-4) continue;
            scale = std::sqrt(scale);

            const auto& offset = evec.col(i) * scale;
            const auto& p1 = mean - offset;
            const auto& p2 = mean + offset;
            marker.points.push_back(to_point(p1));
            marker.points.push_back(to_point(p2));

            marker.color.r = color.r;
            marker.color.g = color.g;
            marker.color.b = color.b;
            dst->markers.push_back(marker);
        }
    };

    const auto& storages = src->getStorages();
    for (const auto& storage : storages) {
        storage->traverse(process_item);
    }
}
}
}

#endif // CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
