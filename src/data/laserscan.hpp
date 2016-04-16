#ifndef LASERSCAN_HPP
#define LASERSCAN_HPP

#include "pointcloud.hpp"
#include <eigen3/Eigen/Core>

namespace ndt {
namespace data {
struct LaserScan : Pointcloud<2> {
    typedef std::shared_ptr<LaserScan> Ptr;

    LaserScan() :
        BaseClass(),
        ranges(nullptr),
        angles(nullptr)
    {
    }

    LaserScan(const LaserScan &other) :
        BaseClass(other),
        ranges(new float[size]),
        angles(new float[size])
    {
        std::memcpy(ranges, other.ranges, sizeof(float) * size);
        std::memcpy(angles, other.angles, sizeof(float) * size);
    }

    LaserScan & operator = (const LaserScan &other)
    {
        if(this != &other) {
            std::size_t former_size = size;
            Pointcloud::operator =(other);
            if(size != former_size) {
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

    inline void resize(const std::size_t _size) override
    {
        if(size != _size) {
            clear();
            BaseClass::resize(_size);
            ranges = new float[_size];
            angles = new float[_size];
        }
        memset(ranges, 0, size * sizeof(float));
        memset(angles,0, size * sizeof(float));
    }

    inline void clear() override
    {
        BaseClass::clear();
        delete[] ranges;
        ranges = nullptr;
        delete[] angles;
        angles = nullptr;
    }

    float     *ranges;
    float     *angles;
};
}
}

#endif // LASERSCAN_HPP
