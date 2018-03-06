#ifndef CSLIBS_NDT_2D_MAP_LOADER_H
#define CSLIBS_NDT_2D_MAP_LOADER_H

#include <cslibs_ndt_2d/serialization/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_2d/serialization/dynamic_maps/occupancy_gridmap.hpp>

#include <nav_msgs/OccupancyGrid.h>

#include <ros/ros.h>
#include <ros/service.h>
#include <std_srvs/Empty.h>

namespace cslibs_ndt_2d {
class NDTMapLoader
{
public:
    NDTMapLoader();
    virtual ~NDTMapLoader();

    void run();

private:
    ros::NodeHandle     nh_;
    ros::Publisher      pub_ndt_grid_;
    ros::Publisher      pub_occ_ndt_grid_;
    ros::ServiceServer  service_;

    cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr          map_ndt_;
    cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr map_occ_ndt_;

    nav_msgs::OccupancyGrid::Ptr                       map_ndt_grid_;
    nav_msgs::OccupancyGrid::Ptr                       map_occ_ndt_grid_;


    bool setup();

    bool resend(std_srvs::Empty::Request &req,
                std_srvs::Empty::Response &res);
};
}

#endif // CSLIBS_NDT_2D_MAP_LOADER_H
