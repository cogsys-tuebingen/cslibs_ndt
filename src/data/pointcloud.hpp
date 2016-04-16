#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <memory>
#include <cstring>

namespace ndt {
namespace data {
template<typename PointT>
struct Pointcloud {
    typedef std::shared_ptr<Pointcloud<PointT>> Ptr;
    typedef PointT             PointType;
    typedef Pointcloud<PointT> BaseClass;

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

    Pointcloud(const Pointcloud &other) :
        size(other.size),
        points(new PointType[size]),
        mask(new char[size])
    {
        std::memcpy(points, other.points, sizeof(PointType) * size);
        std::memcpy(mask, other.mask, size);
    }

    Pointcloud & operator = (const Pointcloud &other)
    {
        if(this != &other) {
            std::size_t former_size = size;
            size = other.size;
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

    virtual void resize(const std::size_t _size)
    {
        clear();
        size = _size;
        points = new PointT[size];
        mask = new char[size];
        memset(points, 0, size * sizeof(PointT));
        memset(mask, INVALID, size);
    }

    virtual void clear()
    {
        size = 0;
        delete[] points;
        points = nullptr;
        delete[] mask;
        mask = nullptr;
    }

    std::size_t size;
    PointT     *points;
    char       *mask;
};
}
}
#endif // POINTCLOUD_HPP
