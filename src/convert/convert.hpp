#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <sensor_msgs/LaserScan.h>

#include "../data/laserscan.hpp"

namespace ndt {
inline void convert(const sensor_msgs::LaserScanConstPtr &msg,
                    ndt::data::LaserScan &scan)
{
    scan.resize(msg->ranges.size());
    std::memcpy(scan.ranges, msg->ranges.data(), sizeof(float) * scan.size);

    float angle = msg->angle_min;
    float angle_incr = msg->angle_increment;
    float range_min = msg->range_min;
    float range_max = mag->range_max;

    for(std::size_t i = 0 ; i < scan.size ; ++i, angle+=angle_incr) {
        float r = scan.ranges[i];
        if(r >= range_min &&
           r <= range_max) {
            Point &p = scan.points[i];
            sincos(angle, &(p(0)), &(p(1)));
            p *= r;
            scan.mask = LaserScan::VALID;
        }
        scan.angles[i] = angle;
    }

    scan.min = Eigen::Vector2d(-range_max, -range_max);
    scan.max = Eigen::Vector2d(range_max, range_max);
}
}

#endif // CONVERT_HPP
