#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <sensor_msgs/LaserScan.h>

#include <ndt/data/laserscan.hpp>

namespace ndt {
namespace conversion {
inline void convert(const sensor_msgs::LaserScanConstPtr &msg,
                    ndt::data::LaserScan &scan,
                    bool inverted = false)
{
    scan.resize(msg->ranges.size());

    float angle = msg->angle_min;
    float angle_incr = msg->angle_increment;
    float range_min = msg->range_min;
    float range_max = msg->range_max;

    for(std::size_t i = 0 ; i < scan.size ; ++i, angle+=angle_incr) {
        std::size_t index;
        if(inverted)
            index = scan.size - 1 - i;
        else
            index = i;

        float r = msg->ranges[index];
        if(r >= range_min &&
           r <= range_max) {
            data::LaserScan::PointType &p = scan.points[index];
            sincos(angle, &(p(1)), &(p(0)));
            p *= r;
            scan.mask[index] = data::LaserScan::VALID;
        }
        scan.ranges[index] = r;
        scan.angles[index] = angle;
    }

    scan.min = Eigen::Vector2d(-range_max, -range_max);
    scan.max = Eigen::Vector2d(range_max, range_max);
}
}
}

#endif // CONVERT_HPP
