#ifndef NDT_GRID_HPP
#define NDT_GRID_HPP

#include <ndt/math/distribution.hpp>

#include <array>
#include <iostream>

namespace ndt {
namespace grid {
template<std::size_t Dim>
class Grid {
public:
    typedef std::shared_ptr<Grid<Dim>>    Ptr;

    typedef std::array<std::size_t, Dim>      SizeType;
    typedef std::array<std::size_t, Dim>      IndexType;
    typedef std::array<double, Dim>           ResolutionType;
    typedef Eigen::Matrix<double, Dim, 1>     PointType;
    typedef math::Distribution<Dim, true>     DistributionType;
    typedef typename DistributionType::MatrixType CovarianceMatrixType;

    Grid() :
        data_size(0),
        data(nullptr)
    {
    }

    Grid(const SizeType       &_size,
            const ResolutionType &_resolution,
            const PointType      &_origin = PointType::Zero()) :
        size(_size),
        resolution(_resolution),
        origin(_origin),
        data_size(1),
        data(nullptr)
    {
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            data_size *= _size[i];
        }

        steps[0] = 1;
        if(Dim > 1) {
            std::size_t max_idx = Dim - 1;
            for(std::size_t i = 1 ; i <= max_idx ; ++i) {
                steps[i] = steps[i-1] * size[i];
            }
        }
        data = new DistributionType[data_size];
    }

    Grid(const Grid &other) :
        size(other.size),
        steps(other.steps),
        resolution(other.resolution),
        origin(other.origin),
        data_size(other.data_size),
        data(new DistributionType[data_size])
    {
        std::memcpy(data, other.data, sizeof(DistributionType) * data_size);
    }

    Grid & operator = (const Grid &other)
    {
        if(this != &other) {
            size = other.size;
            steps = other.steps;
            resolution = other.resolution;
            origin = other.origin;
            std::size_t former_size = data_size;
            data_size = other.data_size;
            if(data_size != former_size) {
                delete [] data;
                data = new DistributionType[data_size];
            }

            std::memcpy(data, other.data, sizeof(DistributionType) * data_size);
        }

        return *this;
    }

    virtual ~Grid()
    {
        delete[] data;
        data = nullptr;
    }


    /// ---------------- META INFORMATION ---------------- ///
    inline SizeType getSize() const
    {
        return size;
    }

    inline void getSize(SizeType &_size) const
    {
        _size = size;
    }

    inline PointType getOrigin() const
    {
        return origin;
    }

    inline ResolutionType getResolution() const
    {
        return resolution;
    }

    inline void getResolution(ResolutionType &_resolution) const
    {
        _resolution = resolution;
    }

    inline bool getIndex(const PointType &_p,
                         IndexType &index)
    {
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            int id = (_p(i) - origin(i)) / resolution[i];
            if(id < 0 || id >= size[i])
                return false;
            index[i] = id;
        }
        return true;
    }

    inline bool checkIndex(const IndexType &_index)
    {
        std::size_t p = pos(_index);
        return p < data_size;
    }


    /// ---------------- DATA ---------------------------- ///
    inline bool add(const PointType &_p)
    {
        int p = pos(_p);
        if(p >= data_size || p < 0)
            return false;
        data[p].add(_p);
        return true;
    }

    inline DistributionType* get(const PointType &_p)
    {
        int p = pos(_p);
        if(p >= data_size || p < 0)
            return nullptr;
        return &data[p];
    }

    inline DistributionType const * get(const PointType &_p) const
    {
        int p = pos(_p);
        if(p >= data_size || p < 0)
            return nullptr;
        return &data[p];
    }

    inline DistributionType const & at(const IndexType &_index) const
    {
        int p = pos(_index);
        if(p >= data_size || p < 0)
            throw std::runtime_error("Out of bounds!");

        return data[p];
    }

    inline DistributionType & at(const IndexType &_index)
    {
        int p = pos(_index);
        if(p >= data_size || p < 0)
            throw std::runtime_error("Out of bounds!");
        return data[p];
    }

private:
    SizeType             size;
    SizeType             steps;
    ResolutionType       resolution;
    PointType            origin;
    std::size_t      data_size;
    DistributionType    *data;

    inline std::size_t pos(const IndexType &_index) {
        std::size_t pos;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            pos += _index[i] * steps[i];
        }
        return pos;
    }

    inline int pos(const PointType &_p) {
        int pos = 0;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            pos += floor((_p(i) - origin(i)) / resolution[i])* steps[i];
        }

        return pos;
    }
};
}
}
#endif // NDT_GRID_HPP
