#ifndef NDT_GRID_HPP
#define NDT_GRID_HPP
#include "../math/distribution.hpp"
#include <array>

namespace ndt {
template<std::size_t Dim>
class NDTGrid {
public:
    typedef std::shared_ptr<NDTGrid<Dim>> Ptr;
    typedef std::array<std::size_t, Dim>  Size;
    typedef std::array<std::size_t, Dim>  Index;
    typedef std::array<double, Dim>       Resolution;
    typedef Eigen::Matrix<double, Dim, 1> Point;
    typedef math::Distribution<Dim, true> Distribution;
    typedef typename Distribution::Matrix Matrix;

    NDTGrid() :
        data_size(0),
        data(nullptr)
    {
    }

    NDTGrid(const Size       &_size,
            const Resolution &_resolution,
            const Point      &_origin = Point::Zero()) :
        size(_size),
        resolution(_resolution),
        origin(_origin),
        data_size(1),
        data(nullptr)
    {
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            data_size *= _size[i];
        }

        steps[Dim - 1] = 1;
        if(Dim > 1) {
            std::size_t max_idx = Dim - 1;
            for(std::size_t i = max_idx ; i > 0 ; --i) {
                steps[i - 1] = steps[i] * size[i];
            }
        }
        data = new Distribution[data_size];
    }

    virtual ~NDTGrid()
    {
        if(data) {
            delete[] data;
        }
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

    inline bool getIndex(const Point &_p,
                         Index &index)
    {
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            int id = (_p(i) - origin(i)) / resolution(i);
            if(id < 0 || id >= size[i])
                return false;
            index[i] = id;
        }
        return true;
    }

    inline bool checkIndex(const Index &_index)
    {
        std::size_t p = pos(_index);
        return p < data_size;
    }


    /// ---------------- DATA ---------------------------- ///
    inline bool add(const Point &_p)
    {
        std::size_t p = pos(_p);
        if(p >= data_size)
            return false;
        data[p].add(_p);
        return true;
    }

    inline double sample(const Point &_p)
    {
        std::size_t p = pos(_p);
        if(p >= data_size)
            return 0.0;
        return data[p].evaluate(_p);
    }

    inline double sample(const Point &_p,
                         Point       &_mean,
                         Matrix      &_inverse_covariance)
    {
        std::size_t p = pos(_p);
        if(p >= data_size)
            return 0.0;
        data[p].getMean(_mean);
        data[p].getInversecovariance(_inverse_covariance);
        return data[p].evaluate(_p);
    }

    inline double sampleNonNormalized(const Point &_p)
    {
        std::size_t p = pos(_p);
        if(p >= data_size)
            return 0.0;
        return data[p].evaluateNonNoramlized(_p);
    }

    inline double sampleNonNormalized(const Point &_p,
                                      Point       &_mean,
                                      Matrix      &_inverse_covariance)
    {
        std::size_t p = pos(_p);
        if(p >= data_size)
            return 0.0;
        data[p].getMean(_mean);
        data[p].getInversecovariance(_inverse_covariance);
        return data[p].evaluateNonNoramlized(_p);
    }

    inline Distribution const & at(const Index &_index) const
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");

        return data[p];
    }

    inline Distribution & at(const Index &_index)
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");
        return data[p];
    }

private:
    Size             size;
    Size             steps;
    Resolution       resolution;
    Point            origin;
    std::size_t      data_size;
    Distribution    *data;

    inline std::size_t pos(const Index &_index) {
        std::size_t pos;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            pos += _index[i] * steps[i];
        }
        return pos;
    }

    inline std::size_t pos(const Point &_p) {
        std::size_t pos = 0;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            pos += (_p(i) - origin(i)) / resolution[i] * steps[i];
        }
        return pos;
    }
};
}
#endif // NDT_GRID_HPP
