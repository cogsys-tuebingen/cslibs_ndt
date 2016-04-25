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
        points_data(size, PointType::Zero()),
        points(points_data.data()),
        mask_data(size, INVALID),
        mask(mask_data.data()),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
    }

    Pointcloud(const std::vector<PointType> &_points) :
        size(_points.size()),
        points_data(_points),
        points(points_data.data()),
        mask_data(size, VALID),
        mask(mask_data.data()),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
        autoAdjustLimits();
    }

    Pointcloud(const Pointcloud &other) :
        size(other.size),
        points_data(other.points_data),
        points(points_data.data()),
        mask_data(other.mask_data),
        mask(mask_data.data()),
        min(other.min),
        max(other.max)
    {
    }

    Pointcloud & operator = (const Pointcloud &other)
    {
        if(this != &other) {
            size = other.size;
            min = other.min;
            max = other.max;
            points_data = other.points_data;
            points = points_data.data();
            mask_data = other.mask_data;
            mask = mask_data.data();
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
        PointType _max = PointType::Constant(std::numeric_limits<double>::lowest());

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
        size = _size;
        points_data.resize(size, PointType::Zero());
        points = points_data.data();
        mask_data.resize(size, INVALID);
        mask = mask_data.data();
    }

    inline virtual void clear()
    {
        size = 0;
        points_data.clear();
        points = nullptr;
        mask_data.clear();
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
    std::vector<PointType> points_data;
    PointType  *points;
    std::vector<char>      mask_data;
    char       *mask;
    PointType   min;
    PointType   max;
};
}
}
#endif // POINTCLOUD_HPP
