#ifndef LASERSCAN_HPP
#define LASERSCAN_HPP

#include <ndt/data/pointcloud.hpp>
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
        ranges_data(other.ranges_data),
        ranges(ranges_data.data()),
        angles_data(other.angles_data),
        angles(angles_data.data())
    {
    }

    LaserScan & operator = (const LaserScan &other)
    {
        if(this != &other) {
            BaseClass::operator =(other);
            ranges_data = other.ranges_data;
            ranges = ranges_data.data();
            angles_data = other.angles_data;
            angles = angles_data.data();
        }
        return *this;
    }


    virtual ~LaserScan()
    {
        clear();
    }

    inline void resize(const std::size_t _size) override
    {
        BaseClass::resize(_size);
        ranges_data.resize(_size, 0.f);
        ranges = ranges_data.data();
        angles_data.resize(_size, 0.f);
        angles = angles_data.data();
    }

    inline void clear() override
    {
        BaseClass::clear();
        ranges_data.clear();
        ranges = nullptr;
        angles_data.clear();
        angles = nullptr;
    }

    std::vector<float> ranges_data;
    float     *ranges;
    std::vector<float> angles_data;
    float     *angles;
};
}
}

#endif // LASERSCAN_HPP
