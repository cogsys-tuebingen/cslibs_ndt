#ifndef CSLIBS_NDT_2D_MAP_LOADER_H
#define CSLIBS_NDT_2D_MAP_LOADER_H

#include <cslibs_ndt_2d/serialization/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/serialization/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_ndt_2d/conversion/probability_gridmap.hpp>
#include <cslibs_ndt_2d/conversion/distributions.hpp>

#include <cslibs_gridmaps/static_maps/algorithms/normalize.hpp>
#include <cslibs_gridmaps/static_maps/conversion/convert_probability_gridmap.hpp>

#include <nav_msgs/OccupancyGrid.h>
#include <visualization_msgs/MarkerArray.h>

#include <ros/ros.h>
#include <ros/service.h>
#include <std_srvs/Empty.h>

namespace cslibs_ndt_2d {
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
    ros::Publisher      pub_ndt_grid_;
    ros::Publisher      pub_occ_ndt_grid_;
    ros::Publisher      pub_ndt_distributions_;
    ros::Publisher      pub_occ_ndt_distributions_;
    ros::ServiceServer  service_;

    nav_msgs::OccupancyGrid::Ptr map_ndt_grid_;
    nav_msgs::OccupancyGrid::Ptr map_occ_ndt_grid_;

    visualization_msgs::MarkerArray::Ptr map_ndt_distributions_;
    visualization_msgs::MarkerArray::Ptr map_occ_ndt_distributions_;

    template <typename T>
    bool setup()
    {
        const std::string path_ndt     = nh_.param<std::string>("path_ndt",     "");
        const std::string path_occ_ndt = nh_.param<std::string>("path_occ_ndt", "");

        const std::string topic_ndt_grid     = nh_.param<std::string>("topic_ndt_grid",     "/map/2d/ndt");
        const std::string topic_occ_ndt_grid = nh_.param<std::string>("topic_occ_ndt_grid", "/map/2d/occ_ndt");

        pub_ndt_grid_     = nh_.advertise<nav_msgs::OccupancyGrid>(topic_ndt_grid,     1);
        pub_occ_ndt_grid_ = nh_.advertise<nav_msgs::OccupancyGrid>(topic_occ_ndt_grid, 1);

        const double ndt_sampling_resolution     = nh_.param<double>("ndt_sampling_resolution",     0.025);
        const double occ_ndt_sampling_resolution = nh_.param<double>("occ_ndt_sampling_resolution", 0.025);

        const std::string topic_ndt_distributions     = nh_.param<std::string>("topic_ndt_distributions",     "/map/2d/ndt/distributions");
        const std::string topic_occ_ndt_distributions = nh_.param<std::string>("topic_occ_ndt_distributions", "/map/2d/occ_ndt/distributions");

        pub_ndt_distributions_     = nh_.advertise<visualization_msgs::MarkerArray>(topic_ndt_distributions,     1);
        pub_occ_ndt_distributions_ = nh_.advertise<visualization_msgs::MarkerArray>(topic_occ_ndt_distributions, 1);

        const ros::Time time    = ros::Time::now();
        const std::string frame = nh_.param<std::string>("map_frame", "/map");

        if (path_ndt != "") {
            typename cslibs_ndt_2d::dynamic_maps::Gridmap<T>::Ptr map_ndt;
            if (!cslibs_ndt_2d::dynamic_maps::loadBinary<T>(path_ndt, map_ndt)) {
                std::cerr << "Could not load ndt 2d map '" << path_ndt << "'." << std::endl;
                return false;
            }

            typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr tmp;
            cslibs_ndt_2d::conversion::from<T>(map_ndt, tmp, ndt_sampling_resolution);

            if (tmp) {
                cslibs_gridmaps::static_maps::algorithms::normalize<T>(*tmp);
                cslibs_gridmaps::static_maps::conversion::from(*tmp, map_ndt_grid_);

                if (map_ndt_grid_) {
                    map_ndt_grid_->header.stamp    = time;
                    map_ndt_grid_->header.frame_id = frame;
                }
            }            

            cslibs_ndt_2d::conversion::from<T>(map_ndt, map_ndt_distributions_, time, frame);
        }
        if (path_occ_ndt != "") {
            typename cslibs_ndt_2d::dynamic_maps::OccupancyGridmap<T>::Ptr map_occ_ndt;
            if (!cslibs_ndt_2d::dynamic_maps::loadBinary<T>(path_occ_ndt, map_occ_ndt)) {
                std::cerr << "Could not load occupancy ndt 2d map '" << path_occ_ndt << "'." << std::endl;
                return false;
            }

            typename cslibs_gridmaps::utility::InverseModel<T>::Ptr ivm(new cslibs_gridmaps::utility::InverseModel<T>(0.5, 0.45, 0.65));
            typename cslibs_gridmaps::static_maps::ProbabilityGridmap<T,T>::Ptr tmp;
            cslibs_ndt_2d::conversion::from<T>(map_occ_ndt, tmp, occ_ndt_sampling_resolution, ivm);

            if (tmp) {
                cslibs_gridmaps::static_maps::algorithms::normalize<T>(*tmp);
                cslibs_gridmaps::static_maps::conversion::from(*tmp, map_occ_ndt_grid_);

                if (map_occ_ndt_grid_) {
                    map_occ_ndt_grid_->header.stamp    = time;
                    map_occ_ndt_grid_->header.frame_id = frame;
                }
            }

            cslibs_ndt_2d::conversion::from<T>(map_occ_ndt, map_occ_ndt_distributions_, ivm, time, frame);
       }

        service_ = nh_.advertiseService(nh_.getNamespace() + "/resend", &NDTMapLoader::resend, this);

        return true;
    }

    bool resend(std_srvs::Empty::Request &req,
                std_srvs::Empty::Response &res)
    {
        if (!map_ndt_grid_ && !map_occ_ndt_grid_) {
            ROS_ERROR_STREAM("What can I say, I have nothing to offer!");
            return false;
        }
        const ros::Time time = ros::Time::now();

        if (map_ndt_grid_) {
            map_ndt_grid_->header.stamp = time;
            pub_ndt_grid_.publish(map_ndt_grid_);
        }
        if (map_occ_ndt_grid_) {
            map_occ_ndt_grid_->header.stamp = time;
            pub_occ_ndt_grid_.publish(map_occ_ndt_grid_);
        }

        if (map_ndt_distributions_) {
            for (auto& m : map_ndt_distributions_->markers)
                m.header.stamp = time;
            pub_ndt_distributions_.publish(map_ndt_distributions_);
        }
        if (map_occ_ndt_distributions_) {
            for (auto& m : map_occ_ndt_distributions_->markers)
                m.header.stamp = time;
            pub_occ_ndt_distributions_.publish(map_occ_ndt_distributions_);
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_MAP_LOADER_H
