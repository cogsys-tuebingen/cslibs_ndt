#include "ndt_map_loader.h"

#include <cslibs_ndt_3d/conversion/pointcloud.hpp>
#include <cslibs_ndt_3d/conversion/distributions.hpp>

namespace cslibs_ndt_3d {
NDTMapLoader::NDTMapLoader() :
    nh_("~")
{
}

NDTMapLoader::~NDTMapLoader()
{
}

void NDTMapLoader::run()
{
    std::cerr << "Starting to set up." << std::endl;
    if (!setup()) {
        std::cerr << "Cannot setup the map loader node." << std::endl;
        ros::shutdown();
    }
    std::cerr << "Setup succesful." << std::endl;

    ros::spin();
}

bool NDTMapLoader::setup()
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
        if (!cslibs_ndt_3d::dynamic_maps::loadBinary(path_ndt, map_ndt_)) {
            std::cerr << "Could not load ndt 3d map '" << path_ndt << "'." << std::endl;
            return false;
        }

        pcl::PointCloud<pcl::PointXYZI>::Ptr points;
        cslibs_ndt_3d::conversion::from(map_ndt_, points);
        map_ndt_means_.reset(new sensor_msgs::PointCloud2);
        pcl::toROSMsg(*points, *map_ndt_means_);
        map_ndt_means_->header.frame_id = "/map";

        cslibs_ndt_3d::conversion::from(map_ndt_, map_ndt_distributions_);
        map_ndt_distributions_->header.frame_id = "/map";
    }
    if (path_occ_ndt != "") {
        if (!cslibs_ndt_3d::dynamic_maps::loadBinary(path_occ_ndt, map_occ_ndt_)) {
            std::cerr << "Could not load occupancy ndt 3d map '" << path_occ_ndt << "'." << std::endl;
            return false;
        }
        cslibs_gridmaps::utility::InverseModel::Ptr ivm(new cslibs_gridmaps::utility::InverseModel(0.5, 0.45, 0.65));
        pcl::PointCloud<pcl::PointXYZI>::Ptr points;
        cslibs_ndt_3d::conversion::from(map_occ_ndt_, points, ivm);
        map_occ_ndt_means_.reset(new sensor_msgs::PointCloud2);
        pcl::toROSMsg(*points, *map_occ_ndt_means_);
        map_occ_ndt_means_->header.frame_id = "/map";

        cslibs_ndt_3d::conversion::from(map_occ_ndt_, map_occ_ndt_distributions_, ivm);
        map_occ_ndt_distributions_->header.frame_id = "/map";
    }

    service_ = nh_.advertiseService(nh_.getNamespace() + "/resend", &NDTMapLoader::resend, this);

    return true;
}

bool NDTMapLoader::resend(std_srvs::Empty::Request  &req,
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
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "NDTMapLoader");

    cslibs_ndt_3d::NDTMapLoader ml;
    ml.run();

    return 0;
}
