#ifndef MULTI_GRID_HPP
#define MULTI_GRID_HPP

#include "grid.hpp"
#include "mask.hpp"
#include "../data/pointcloud.hpp"

namespace ndt {
template<std::size_t Dim>
class NDTMultiGrid {
public:
    typedef std::shared_ptr<NDTMultiGrid<Dim>>  Ptr;
    typedef typename NDTGrid<Dim>::Size         Size;
    typedef typename NDTGrid<Dim>::Index        Index;
    typedef typename NDTGrid<Dim>::Resolution   Resolution;
    typedef typename NDTGrid<Dim>::Point        Point;
    typedef typename NDTGrid<Dim>::Distribution Distribution;
    typedef typename NDTGrid<Dim>::Matrix       Matrix;


    NDTMultiGrid() :
        data_size(0),
        data(nullptr)
    {
    }
    /// NOTICE :
    /// Displacement of 0.25 * resolution results in 0.5 displacement from grids to each other
    /// Origin of Multigrid is in the center of the combined grid area.

    /// maybe point cloud constructor

    NDTMultiGrid(const Size       &_size,
                 const Resolution &_resolution,
                 const Point      &_origin = Point::Zero()) :
        size(_size),
        resolution(_resolution),
        origin(_origin),
        data_size(mask.rows),
        data(new NDTGrid<Dim>[data_size])
    {
        Resolution offsets;
        for(std::size_t i = 0 ; i < Dim; ++i) {
            offsets[i] = +_resolution[i] * 0.25;
        }

        for(std::size_t i = 0 ; i < data_size ; ++i) {
            Point o = origin;
            for(std::size_t j = 0 ; j < Dim ; ++j) {
                o(j) += mask[i * mask.cols + j] * offsets[j];
            }
            data[i] = NDTGrid<Dim>(_size, _resolution, o);
        }

        steps[0] = 1;
        if(Dim > 1) {
            std::size_t max_idx = Dim - 1;
            for(std::size_t i = 1 ; i <= max_idx ; ++i) {
                steps[i] = steps[i-1] * 2;
            }
        }
    }

    NDTMultiGrid(const NDTMultiGrid &other) :
        size(other.size),
        resolution(other.resolution),
        origin(other.origin),
        data_size(other.data_size),
        data(new NDTGrid<Dim>[data_size])
    {
        std::memcpy(data, other.data, sizeof(NDTGrid<Dim>) * data_size);
    }

    NDTMultiGrid & operator = (const NDTMultiGrid &other)
    {
        if(this != &other) {
            std::size_t former_size = size;
            size = other.size;
            resolution = other.resolution;
            origin = other.origin;
            data_size = other.data_size;
            if(size != former_size) {
                delete [] data;
                data = new NDTGrid<Dim>[data_size];
            }
            std::memcpy(data, other.data, sizeof(NDTGrid<Dim>) * data_size);
        }
        return *this;
    }

    virtual ~NDTMultiGrid()
    {
        delete[] data;
        data = nullptr;
    }

    /// ---------------- META INFORMATION ---------------- ///
    inline Size getSize() const
    {
        return size;
    }

    inline void getSize(Size &_size) const
    {
        _size = size;
    }

    inline Resolution getResolution() const
    {
        return resolution;
    }

    inline void getResolution(Resolution &_resolution) const
    {
        _resolution = resolution;
    }

    inline bool checkIndex(const Index &_index)
    {
        std::size_t p = pos(_index);
        return p < data_size;
    }

    inline Point getOrigin() const
    {
        return origin;
    }

    /// ---------------- DATA ---------------------------- ///
    inline bool add(const Point &_p)
    {
        bool result = false;
        for(std::size_t i = 0 ; i < data_size; ++i) {
            result |= data[i].add(_p);
        }
        return result;
    }

    inline bool add(const data::Pointcloud<Dim> &_p)
    {
        bool result = false;
        for(std::size_t i = 0 ; i < _p.size ; ++i) {
            if(_p.mask[i] == data::Pointcloud<Dim>::VALID) {
                for(std::size_t j = 0 ; j < data_size; ++j)
                    result |= data[j].add(_p.points[i]);
            }
        }
        return result;
    }

    inline double sample(const Point &_p)
    {
        double result = 0.0;
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sample(_p);
        }
        return result;
    }

    inline double sample(const Point &_p,
                         Point       &_mean,
                         Matrix      &_inverse_covariance)
    {
        _mean = Point::Zero();
        _inverse_covariance = Matrix::Zero();

        double result = 0.0;
        Point   mean = Point::Zero();
        Matrix  inverse_covariance = Matrix::Zero();
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sample(_p, mean, inverse_covariance);
            _mean  += mean;
            _inverse_covariance += inverse_covariance;
        }
        return result;
    }

    inline double sample(const Point &_p,
                         Point       &_mean,
                         Matrix      &_inverse_covariance,
                         Point       &_q)
    {
        _mean = Point::Zero();
        _inverse_covariance = Matrix::Zero();

        double result = 0.0;
        Point   mean = Point::Zero();
        Point   q = Point::Zero();
        Matrix  inverse_covariance = Matrix::Zero();
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sample(_p, mean, inverse_covariance, q);
            _mean  += mean;
            _q     += q;
            _inverse_covariance += inverse_covariance;
        }
        return result;
    }

    inline double sampleNonNormalized(const Point &_p)
    {
        double result = 0.0;
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sampleNonNormalized(_p);
        }
        return result;
    }

    inline double sampleNonNormalized(const Point &_p,
                                      Point       &_mean,
                                      Matrix      &_inverse_covariance,
                                      Point       &_q)
    {
        _mean = Point::Zero();
        _inverse_covariance = Matrix::Zero();

        double result = 0.0;
        Point   mean = Point::Zero();
        Matrix  inverse_covariance = Matrix::Zero();
        Point   q = Point::Zero();
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sampleNonNormalized(_p, mean, inverse_covariance, q);
            _mean  += mean;
            _q     += q;
            _inverse_covariance += inverse_covariance;
        }
        return result;
    }

    inline NDTGrid<Dim> const & at(const Index &_index) const
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");

        return data[p];
    }

    inline NDTGrid<Dim> & at(const Index &_index)
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");

        return data[p];
    }

private:
    Size          size;
    Size          steps;
    Resolution    resolution;
    Point         origin;
    std::size_t   data_size;
    NDTGrid<Dim> *data;

    Mask<Dim>     mask;

    inline std::size_t pos(const Index &_index) {
        std::size_t p = 0;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            p += steps[i] * _index[i];
        }
        return p;
    }

};
}
#endif // MULTI_GRID_HPP
