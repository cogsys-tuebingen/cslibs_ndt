#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <sensor_msgs/LaserScan.h>
#include <ndt/data/laserscan.hpp>

namespace ndt {
namespace conversion {
inline void convert(const sensor_msgs::LaserScanConstPtr &_msg,
                    ndt::data::LaserScan &_scan,
                    const float _range_min = -1.f,
                    const float _range_max = -1.f)
{
    float range_min = _range_min;
    float range_max = _range_max;
    if(range_min < 0.f)
        range_min = _msg->range_min;
    if(range_max < 0.f)
        range_max = _msg->range_max;

    _scan.resize(_msg->ranges.size());
    float angle = _msg->angle_min;
    float angle_incr = _msg->angle_increment;

    for(std::size_t i = 0 ; i < _scan.size ; ++i, angle+=angle_incr) {
        float r = _msg->ranges[i];
        if(r >= range_min &&
           r <= range_max) {
            data::LaserScan::PointType &p = _scan.points[i];
            sincos(angle, &(p(1)), &(p(0)));
            p *= r;
            _scan.mask[i] = data::LaserScan::VALID;
        } else {
            _scan.mask[i] = data::LaserScan::INVALID;
        }
        _scan.ranges[i] = r;
        _scan.angles[i] = angle;
    }

    _scan.min = Eigen::Vector2d(-range_max, -range_max);
    _scan.max = Eigen::Vector2d(range_max, range_max);
}
}
}

#endif // CONVERT_HPP
