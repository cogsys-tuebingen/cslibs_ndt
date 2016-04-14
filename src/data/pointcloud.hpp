#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <memory>
#include <cstring>

namespace ndt {
namespace data {
template<typename PointT>
struct Pointcloud {
    typedef std::shared_ptr<Pointcloud<PointT>> Ptr;
    typedef PointT PointType;
    enum EntyValidity {INVALID = 0, VALID = 1};

    Pointcloud() :
        size(0),
        points(nullptr),
        mask(nullptr)
    {
    }

    Pointcloud(const std::size_t _size) :
        size(_size),
        points(new PointT[size]),
        mask(new char[size])
    {
        memset(points, 0, size * sizeof(PointT));
        memset(mask, INVALID, size);
    }


    virtual ~Pointcloud()
    {
        clear();
    }

    virtual void resize(const std::size_t _size)
    {
        clear();
        size = _size;
        points = new PointT[size];
        mask = new char[size];
    }

    virtual void clear()
    {
        size = 0;
        if(points != nullptr)
            delete[] points;
        if(mask != nullptr)
            delete[] mask;
    }

    std::size_t size;
    PointT     *points;
    char       *mask;
};
}
}
#endif // POINTCLOUD_HPP
