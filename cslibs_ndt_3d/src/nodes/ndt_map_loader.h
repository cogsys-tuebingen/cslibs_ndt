#ifndef CSLIBS_NDT_3D_MAP_LOADER_H
#define CSLIBS_NDT_3D_MAP_LOADER_H

#include <cslibs_ndt_3d/serialization/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/serialization/dynamic_maps/occupancy_gridmap.hpp>

#include <sensor_msgs/PointCloud2.h>
#include <cslibs_ndt_3d/DistributionArray.h>

#include <ros/ros.h>
#include <ros/service.h>
#include <std_srvs/Empty.h>

namespace cslibs_ndt_3d {
class NDTMapLoader
{
public:
    NDTMapLoader();
    virtual ~NDTMapLoader();

    void run();

private:
    ros::NodeHandle     nh_;
    ros::Publisher      pub_ndt_means_;
    ros::Publisher      pub_occ_ndt_means_;
    ros::Publisher      pub_ndt_distributions_;
    ros::Publisher      pub_occ_ndt_distributions_;
    ros::ServiceServer  service_;

    cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr          map_ndt_;
    cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr map_occ_ndt_;

    sensor_msgs::PointCloud2::Ptr                      map_ndt_means_;
    sensor_msgs::PointCloud2::Ptr                      map_occ_ndt_means_;

    cslibs_ndt_3d::DistributionArray::Ptr              map_ndt_distributions_;
    cslibs_ndt_3d::DistributionArray::Ptr              map_occ_ndt_distributions_;

    bool setup();

    bool resend(std_srvs::Empty::Request &req,
                std_srvs::Empty::Response &res);
};
}

#endif // CSLIBS_NDT_3D_MAP_LOADER_H
