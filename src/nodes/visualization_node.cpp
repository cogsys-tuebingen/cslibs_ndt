#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <chrono>

#include <ndt/conversion/convert.hpp>
#include <ndt/data/laserscan.hpp>
#include <ndt/grid/multi_grid.hpp>
#include <ndt/visualization/visualize.hpp>

namespace ndt {
struct ScanVisualizerNode {
    typedef ndt::grid::MultiGrid<2> Grid2DType;

    ros::NodeHandle   nh;
    ros::Subscriber   sub;
    double            resolution;
    cv::Mat           display;

    ScanVisualizerNode() :
        nh("~"),
        resolution(1.0)
    {
        std::string topic = "/scan";
        nh.getParam("resolution", resolution);
        nh.getParam("topic", topic);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic, 1, &ScanVisualizerNode::laserscan, this);

    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();
        data::LaserScan scan;
        conversion::convert(msg, scan);
        /// make a grid
        data::LaserScan::PointType range = scan.range();
        Grid2DType::Size       size = {static_cast<std::size_t>(range(0) / resolution),
                                        static_cast<std::size_t>(range(1) / resolution)};
        Grid2DType::Resolution res = {resolution, resolution};

        ndt::MultiGrid2DType grid(size, res, scan.min);
        std::vector<data::LaserScan::PointType> points;
        for(std::size_t i = 0 ; i < scan.size ; ++i) {
            if(scan.mask[i] == data::LaserScan::VALID) {
                points.push_back(scan.points[i]);
                if(!grid.add(scan.points[i]))
                    std::cerr << "Failed to add point [" << scan.points[i] << "]" << std::endl;
            }
        }

        std::chrono::duration<double> elapsed =
                std::chrono::system_clock::now() - start;
        std::cout << "Took : " <<  elapsed.count() * 1000.0 << "ms" << std::endl;

        /// render the grid
        display = cv::Mat (500,500, CV_8UC3, cv::Scalar());
        ndt::renderNDTGrid(grid, scan.min, scan.max, display);
        ndt::renderPoints(points, grid.getSize(), grid.getResolution(), display);
    }

};
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_visualization_node");
    ndt::ScanVisualizerNode sn;
    while(ros::ok()) {
        if(!sn.display.empty()) {
            cv::imshow("ndt", sn.display);
            cv::waitKey(19);
        }
        ros::Rate(30).sleep();
        ros::spinOnce();
    }

    return 0;
}
