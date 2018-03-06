#include "ndt_map_loader.h"

#include <cslibs_ndt_2d/conversion/probability_gridmap.hpp>
#include <cslibs_gridmaps/static_maps/algorithms/normalize.hpp>
#include <cslibs_gridmaps/static_maps/conversion/convert_probability_gridmap.hpp>

namespace cslibs_ndt_2d {
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

    const std::string topic_ndt_grid     = nh_.param<std::string>("topic_ndt_grid",     "/map/2d/ndt");
    const std::string topic_occ_ndt_grid = nh_.param<std::string>("topic_occ_ndt_grid", "/map/2d/occ_ndt");

    pub_ndt_grid_     = nh_.advertise<nav_msgs::OccupancyGrid>(topic_ndt_grid,     1);
    pub_occ_ndt_grid_ = nh_.advertise<nav_msgs::OccupancyGrid>(topic_occ_ndt_grid, 1);

    const double ndt_sampling_resolution     = nh_.param<double>("ndt_sampling_resolution",     0.025);
    const double occ_ndt_sampling_resolution = nh_.param<double>("occ_ndt_sampling_resolution", 0.025);

    if (path_ndt != "") {
        if (!cslibs_ndt_2d::dynamic_maps::loadBinary(path_ndt, map_ndt_)) {
            std::cerr << "Could not load ndt 2d map '" << path_ndt << "'." << std::endl;
            return false;
        }

        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr tmp;
        cslibs_ndt_2d::conversion::from(map_ndt_, tmp, ndt_sampling_resolution);

        if (tmp) {
            cslibs_gridmaps::static_maps::algorithms::normalize<double>(*tmp);
            cslibs_gridmaps::static_maps::conversion::from(tmp, map_ndt_grid_);

            if (map_ndt_grid_)
                map_ndt_grid_->header.frame_id = "/map";
        }
    }
    if (path_occ_ndt != "") {
        if (!cslibs_ndt_2d::dynamic_maps::loadBinary(path_occ_ndt, map_occ_ndt_)) {
            std::cerr << "Could not load occupancy ndt 2d map '" << path_occ_ndt << "'." << std::endl;
            return false;
        }

        cslibs_gridmaps::utility::InverseModel::Ptr ivm(new cslibs_gridmaps::utility::InverseModel(0.5, 0.45, 0.65));
        cslibs_gridmaps::static_maps::ProbabilityGridmap::Ptr tmp;
        cslibs_ndt_2d::conversion::from(map_occ_ndt_, tmp, occ_ndt_sampling_resolution, ivm);

        if (tmp) {
            cslibs_gridmaps::static_maps::algorithms::normalize<double>(*tmp);
            cslibs_gridmaps::static_maps::conversion::from(tmp, map_occ_ndt_grid_);

            if (map_occ_ndt_grid_)
                map_occ_ndt_grid_->header.frame_id = "/map";
        }
    }

    service_ = nh_.advertiseService(nh_.getNamespace() + "/resend", &NDTMapLoader::resend, this);

    return true;
}

bool NDTMapLoader::resend(std_srvs::Empty::Request &req,
                          std_srvs::Empty::Response &res)
{
    if (!map_ndt_grid_ && !map_occ_ndt_grid_) {
        ROS_ERROR_STREAM("What can I say, I have nothing to offer!");
        return false;
    }

    if (map_ndt_grid_) {
        map_ndt_grid_->header.stamp = ros::Time::now();
        pub_ndt_grid_.publish(map_ndt_grid_);
    }
    if (map_occ_ndt_grid_) {
        map_occ_ndt_grid_->header.stamp = ros::Time::now();
        pub_occ_ndt_grid_.publish(map_occ_ndt_grid_);
    }

    return true;
}
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "NDTMapLoader");

    cslibs_ndt_2d::NDTMapLoader ml;
    ml.run();

    return 0;
}
