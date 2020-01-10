#ifndef CSLIBS_NDT_3D_MAP_LOADER_H
#define CSLIBS_NDT_3D_MAP_LOADER_H

#include <cslibs_ndt_3d/serialization/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/serialization/dynamic_maps/occupancy_gridmap.hpp>

#include <sensor_msgs/PointCloud2.h>
#include <cslibs_ndt_3d/DistributionArray.h>

#include <cslibs_ndt_3d/conversion/sensor_msgs_pointcloud2.hpp>
#include <cslibs_ndt_3d/conversion/distributions.hpp>

#include <ros/ros.h>
#include <ros/service.h>
#include <std_srvs/Empty.h>

namespace cslibs_ndt_3d {
class NDTMapLoader
{
public:
    NDTMapLoader() :
        nh_("~")
    {
    }
    virtual ~NDTMapLoader()
    {
    }

    void run()
    {
        std::cerr << "Starting to set up." << std::endl;
        if (!setup<double>() && !setup<float>()) {
            std::cerr << "Cannot setup the map loader node." << std::endl;
            ros::shutdown();
        }
        std::cerr << "Setup succesful." << std::endl;

        ros::spin();
    }

private:
    ros::NodeHandle     nh_;
    ros::Publisher      pub_ndt_means_;
    ros::Publisher      pub_occ_ndt_means_;
    ros::Publisher      pub_ndt_distributions_;
    ros::Publisher      pub_occ_ndt_distributions_;
    ros::ServiceServer  service_;

    sensor_msgs::PointCloud2::Ptr         map_ndt_means_;
    sensor_msgs::PointCloud2::Ptr         map_occ_ndt_means_;

    cslibs_ndt_3d::DistributionArray::Ptr map_ndt_distributions_;
    cslibs_ndt_3d::DistributionArray::Ptr map_occ_ndt_distributions_;

    template <typename T>
    bool setup()
    {
        const std::string path_ndt     = nh_.param<std::string>("path_ndt",     "");
        const std::string path_occ_ndt = nh_.param<std::string>("path_occ_ndt", "");

        const std::string topic_ndt_means     = nh_.param<std::string>("topic_ndt_means",     "/map/3d/ndt/means");
        const std::string topic_occ_ndt_means = nh_.param<std::string>("topic_occ_ndt_means", "/map/3d/occ_ndt/means");

        pub_ndt_means_     = nh_.advertise<sensor_msgs::PointCloud2>(topic_ndt_means,     1);
        pub_occ_ndt_means_ = nh_.advertise<sensor_msgs::PointCloud2>(topic_occ_ndt_means, 1);

        const std::string topic_ndt_distributions     = nh_.param<std::string>("topic_ndt_distributions",     "/map/3d/ndt/distributions");
        const std::string topic_occ_ndt_distributions = nh_.param<std::string>("topic_occ_ndt_distributions", "/map/3d/occ_ndt/distributions");

        pub_ndt_distributions_     = nh_.advertise<cslibs_ndt_3d::DistributionArray>(topic_ndt_distributions,     1);
        pub_occ_ndt_distributions_ = nh_.advertise<cslibs_ndt_3d::DistributionArray>(topic_occ_ndt_distributions, 1);

        if (path_ndt != "") {
            typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr map_ndt;
            if (!cslibs_ndt_3d::dynamic_maps::loadBinary<T>(path_ndt, map_ndt)) {
                std::cerr << "Could not load ndt 3d map '" << path_ndt << "'." << std::endl;
                return false;
            }

            map_ndt_means_.reset(new sensor_msgs::PointCloud2);
            cslibs_ndt_3d::conversion::from<T>(map_ndt, *map_ndt_means_);
            map_ndt_means_->header.frame_id = "/map";

            //cslibs_ndt_3d::conversion::from<T>(map_ndt, map_ndt_distributions_);
            if (map_ndt_distributions_)
                map_ndt_distributions_->header.frame_id = "/map";
        }
        if (path_occ_ndt != "") {
            typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr map_occ_ndt;
            if (!cslibs_ndt_3d::dynamic_maps::loadBinary<T>(path_occ_ndt, map_occ_ndt)) {
                std::cerr << "Could not load occupancy ndt 3d map '" << path_occ_ndt << "'." << std::endl;
                return false;
            }

            typename cslibs_gridmaps::utility::InverseModel<T>::Ptr ivm(new cslibs_gridmaps::utility::InverseModel<T>(0.5, 0.45, 0.65));
            map_occ_ndt_means_.reset(new sensor_msgs::PointCloud2);
            cslibs_ndt_3d::conversion::from<T>(map_occ_ndt, *map_occ_ndt_means_, ivm);
            map_occ_ndt_means_->header.frame_id = "/map";

            //cslibs_ndt_3d::conversion::from<T>(map_occ_ndt, map_occ_ndt_distributions_, ivm);
            if (map_occ_ndt_distributions_)
                map_occ_ndt_distributions_->header.frame_id = "/map";
        }

        service_ = nh_.advertiseService(nh_.getNamespace() + "/resend", &NDTMapLoader::resend, this);

        return true;
    }

    bool resend(std_srvs::Empty::Request &req,
                std_srvs::Empty::Response &res)
    {
        if (!map_ndt_means_ && !map_occ_ndt_means_ && !map_ndt_distributions_ && !map_occ_ndt_distributions_) {
            ROS_ERROR_STREAM("What can I say, I have nothing to offer!");
            return false;
        }

        if (map_ndt_means_) {
            map_ndt_means_->header.stamp = ros::Time::now();
            pub_ndt_means_.publish(map_ndt_means_);
        }
        if (map_occ_ndt_means_) {
            map_occ_ndt_means_->header.stamp = ros::Time::now();
            pub_occ_ndt_means_.publish(map_occ_ndt_means_);
        }

        if (map_ndt_distributions_) {
            map_ndt_distributions_->header.stamp = ros::Time::now();
            pub_ndt_distributions_.publish(map_ndt_distributions_);
        }
        if (map_occ_ndt_distributions_) {
            map_occ_ndt_distributions_->header.stamp = ros::Time::now();
            pub_occ_ndt_distributions_.publish(map_occ_ndt_distributions_);
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_3D_MAP_LOADER_H
