#ifndef LASERSCAN_HPP
#define LASERSCAN_HPP

#include "pointcloud.hpp"
#include <eigen3/Eigen/Core>

namespace ndt {
namespace data {
struct LaserScan : Pointcloud<Eigen::Vector2d> {
    typedef std::shared_ptr<LaserScan> Ptr;

    LaserScan() :
        ranges(nullptr),
        angles(nullptr)
    {
    }

    virtual ~LaserScan()
    {
        clear();
    }

    void clear() override
    {
        size = 0;
        if(ranges != nullptr)
            delete[] ranges;
        if(angles != nullptr)
            delete[] angles;
    }

    float      *ranges;
    float      *angles;
};
}
}

#endif // LASERSCAN_HPP
