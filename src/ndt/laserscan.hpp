#ifndef LASERSCAN_HPP
#define LASERSCAN_HPP

#include <memory>
#include <sensor_msgs/LaserScan.h>
#include "types.hpp"

namespace ndt {
struct LaserScan {
    typedef std::shared_ptr<LaserScan> Ptr;

    enum EntyValidity {INVALID = 0, VALID = 1};

    LaserScan(const sensor_msgs::LaserScanConstPtr &msg) :
        size(msg->ranges.size()),
        ranges(new float[size]),
        angles(new float[size]),
        points(new Point[size]),
        mask(new char[size])
    {
        std::memcpy(ranges, msg->ranges.data(), sizeof(float) * size);
        float angle = msg->angle_min;
        float angle_incr = msg->angle_increment;
        float range_min = msg->range_min;
        float range_max = msg->range_max;

        for(std::size_t i = 0 ; i < size; ++i) {
            float range = ranges[i];
            if(range >= range_min &&
                    range <= range_max) {
                Point &p = points[i];
                sincos(angle, &(p(0)), &(p(2)));
                p *= range;
                mask[i] = VALID;
            } else {
                mask[i] = INVALID;
            }
            angle += angle_incr;
        }
    }

    LaserScan() :
        size(0),
        ranges(nullptr),
        angles(nullptr),
        points(nullptr)
    {
    }

    virtual ~LaserScan()
    {
        clear();
    }

    void clear()
    {
        size = 0;
        if(ranges != nullptr)
            delete[] ranges;
        if(angles != nullptr)
            delete[] angles;
        if(points != nullptr)
            delete[] points;
        if(mask != nullptr)
            delete[] mask;
    }

    std::size_t size;
    float      *ranges;
    float      *angles;
    Point      *points;
    char       *mask;
};
}

#endif // LASERSCAN_HPP
