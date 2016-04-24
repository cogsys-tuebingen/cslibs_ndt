#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <memory>
#include <cstring>
#include <eigen3/Eigen/Core>
#include <vector>

namespace ndt {
namespace data {
template<std::size_t Dim>
struct Pointcloud {
    typedef std::shared_ptr<Pointcloud<Dim>> Ptr;

    typedef Eigen::Matrix<double, Dim, 1>    PointType;
    typedef Pointcloud<Dim>                  BaseClass;

    enum EntyValidity {INVALID = 0, VALID = 1};

    Pointcloud() :
        size(0),
        points(nullptr),
        mask(nullptr),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
    }

    Pointcloud(const std::size_t _size) :
        size(_size),
        points(new PointType[size]),
        mask(new char[size]),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
        std::memset(points, 0, size * sizeof(PointType));
        std::memset(mask, INVALID, size);
    }

    Pointcloud(const std::vector<PointType> &_points) :
        size(_points.size()),
        points(new PointType[size]),
        mask(new char[size]),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
        std::memcpy(points, _points.data(), size * sizeof(PointType));
        std::memset(mask, VALID, size);
        autoAdjustLimits();
    }

    Pointcloud(const Pointcloud &other) :
        size(other.size),
        points(new PointType[size]),
        mask(new char[size]),
        min(other.min),
        max(other.max)
    {
        std::memcpy(points, other.points, sizeof(PointType) * size);
        std::memcpy(mask, other.mask, size);
    }

    Pointcloud & operator = (const Pointcloud &other)
    {
        if(this != &other) {
            std::size_t former_size = size;
            size = other.size;
            min = other.min;
            max = other.max;
            if(size != former_size) {
                delete [] points;
                points = new PointType[size];
            }
            if(size != former_size) {
                delete [] mask;
                mask = new char[size];
            }

            std::memcpy(points, other.points, sizeof(PointType) * size);
            std::memcpy(mask, other.mask, size);
        }
        return *this;
    }

    virtual ~Pointcloud()
    {
        clear();
    }

    inline void autoAdjustLimits()
    {
        if(size == 0)
            return;

        PointType _min = PointType::Constant(std::numeric_limits<double>::max());
        PointType _max = PointType::Constant(std::numeric_limits<double>::min());

        for(std::size_t i = 0 ; i < size ; ++i) {
            if(mask[i] == VALID) {
                for(std::size_t j = 0 ; j < Dim ; ++j) {
                    if(points[i](j) < _min(j))
                        _min(j) = points[i](j);
                    if(points[i](j) > _max(j))
                        _max(j) = points[i](j);
                }
            }
        }

        min = _min;
        max = _max;
    }


    inline virtual void resize(const std::size_t _size)
    {
        if(size != _size) {
            clear();
            size = _size;
            points = new PointType[size];
            mask = new char[size];
        }
        memset(points, 0, size * sizeof(PointType));
        memset(mask, INVALID, size);
    }

    inline virtual void clear()
    {
        size = 0;
        delete[] points;
        points = nullptr;
        delete[] mask;
        mask = nullptr;
        min = PointType::Zero();
        max = PointType::Zero();
    }

    inline void range(PointType &_range) const
    {
        _range = max - min;
    }

    inline PointType range() const
    {
        return max - min;
    }

    std::size_t size;
    PointType  *points;
    char       *mask;
    PointType   min;
    PointType   max;
};
}
}
#endif // POINTCLOUD_HPP
