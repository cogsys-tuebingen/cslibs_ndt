#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include "../ndt/laserscan.hpp"
#include "../ndt/ndt_multi_grid.h"
#include <eigen3/Eigen/Jacobi>
#include "../optional/visualize.hpp"

namespace ndt {
struct ScanVisualizerNode {
    ros::NodeHandle   nh;
    ros::Subscriber   sub;

    double            resolution;
    double            size;
    double            margin;

    std::size_t       conv_iter;
    double            conv_eps;

    NDTMultiGrid::Ptr ndt_grid;
    std::vector<Point> points;

    ScanVisualizerNode() :
        nh("~"),
        resolution(1.0),
        size(0.0),
        margin(0.1),
        conv_iter(100),
        conv_eps(1e-3)
    {
        std::string topic("/scan");
        nh.getParam("topic", topic);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic, 1, &ScanVisualizerNode::laserscan, this);

    }

    void laserscan(const sensor_msgs::LaserScan::ConstPtr &msg)
    {
        LaserScan scan(msg);
        size = msg->range_max;
        ndt_grid.reset(new NDTMultiGrid(size + margin,
                                        size + margin,
                                        resolution));
        insert(scan, ndt_grid);


        cv::Mat display;
        ndt::renderNDTGrid(*ndt_grid, ndt::Point(-size, -size), ndt::Point(size, size), 0.01, display, false);
        ndt::renderPoints(points, display);
        cv::imshow("display", display);
        cv::waitKey(19);

    }

    void insert(const LaserScan   &scan,
                NDTMultiGrid::Ptr &ndt_grid)
    {
        points.clear();
        for(std::size_t i = 0 ; i < scan.size ; ++i) {
            if(scan.mask[i] == LaserScan::VALID) {
                ndt_grid->add(scan.points[i]);
                points.push_back(scan.points[i]);
            }
        }
    }
};
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "scan_matcher");
    ndt::ScanVisualizerNode sn;

    ros::spin();

    return 0;
}
