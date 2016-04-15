#ifndef MULTI_GRID_HPP
#define MULTI_GRID_HPP

#include "grid.hpp"
#include "mask.hpp"

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

    NDTMultiGrid(const Size       &_size,
                 const Resolution &_resolution) :
        size(_size),
        resolution(_resolution),
        data_size(mask.rows),
        data(new NDTGrid<Dim>[data_size])
    {
        Resolution offsets;
        for(std::size_t i = 0 ; i < Dim; ++i) {
            offsets[i] = _resolution[i] / 2.0;
        }

        for(std::size_t i = 0 ; i < data_size ; ++i) {
            Point origin;
            for(std::size_t j = 0 ; j < Dim ; ++j) {
                origin(j) = mask[i * mask.cols] * offsets[j];
            }
            data[i] = NDTGrid<Dim>(_size, _resolution, origin);
        }

        steps[Dim - 1] = 1;
        if(Dim > 1) {
            std::size_t max_idx = Dim - 1;
            for(std::size_t i = max_idx ; i > 0 ; --i) {
                steps[i - 1] = steps[i] * 2;
            }
        }
    }

    virtual ~NDTMultiGrid()
    {
        if(data)
            delete[] data;
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

    /// ---------------- DATA ---------------------------- ///
    inline bool add(const Point &_p)
    {
        bool result = false;
        for(std::size_t i = 0 ; i < data_size; ++i) {
            result |= data[i].add(_p);
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

    inline double sampleNonNormalize(const Point &_p)
    {
        double result = 0.0;
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sampleNonNormalized(_p);
        }
        return result;
    }

    inline double sampleNonNormalized(const Point &_p,
                                      Point       &_mean,
                                      Matrix      &_inverse_covariance)
    {
        _mean = Point::Zero();
        _inverse_covariance = Matrix::Zero();

        double result = 0.0;
        Point   mean = Point::Zero();
        Matrix  inverse_covariance = Matrix::Zero();
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            result += data[i].sampleNonNormalized(_p, mean, inverse_covariance);
            _mean  += mean;
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
