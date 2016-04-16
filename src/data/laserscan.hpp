#ifndef LASERSCAN_HPP
#define LASERSCAN_HPP

#include "pointcloud.hpp"
#include <eigen3/Eigen/Core>

namespace ndt {
namespace data {
struct LaserScan : Pointcloud<Eigen::Vector2d> {
    typedef std::shared_ptr<LaserScan> Ptr;


    LaserScan() :
        BaseClass(),
        ranges(nullptr),
        angles(nullptr),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
    }

    LaserScan(const LaserScan &other) :
        BaseClass(other),
        ranges(new float[size]),
        angles(new float[size]),
        min(other.min),
        max(other.max)
    {
        std::memcpy(ranges, other.ranges, sizeof(float) * size);
        std::memcpy(angles, other.angles, sizeof(float) * size);
    }

    LaserScan & operator = (const LaserScan &other)
    {
        if(this != &other) {
            std::size_t former_size = size;
            Pointcloud::operator =(other);
            if(ranges && size != former_size) {
                    delete [] ranges;
                    ranges = new float[size];
            }
            if(angles && size != former_size) {
                delete [] angles;
                angles = new float[size];
            }
            std::memcpy(ranges, other.ranges, sizeof(float) * size);
            std::memcpy(angles, other.angles, sizeof(float) * size);
            min = other.min;
            max = other.max;
        }
        return *this;
    }


    virtual ~LaserScan()
    {
        clear();
    }

    void resize(const std::size_t _size) override
    {
        clear();
        BaseClass::resize(_size);
        ranges = new float[_size];
        angles = new float[_size];
    }

    void clear() override
    {
        BaseClass::clear();
        if(ranges != nullptr)
            delete[] ranges;
        if(angles != nullptr)
            delete[] angles;
    }

    float     *ranges;
    float     *angles;
    PointType  min;
    PointType  max;
};
}
}

#endif // LASERSCAN_HPP
