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
    typedef std::array<double>            Resolution;
    typedef Eigen::Matrix<double, Dim, 1> Point;
    typedef math::Distribution<Dim>       Distribution;
    typedef Distribution::Matrix          Matrix;

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
            int id = (point(i) - origin(i)) / resolution(i);
            if(id < 0)
                return false;
            index[i] = id;
        }
        return true;
    }

    inline bool checkIndex(const Index &_index)
    {
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            if(_index >= size[i])
                return false;
        }
        return true;
    }


    /// ---------------- DATA ---------------------------- ///

    inline bool add(const Point &_p)
    {

    }

    inline double sample(const Point &_p)
    {

    }

    inline double sample(const Point &_p,
                         Point       &_mean,
                         Matrix      &_covariance)
    {

    }

    inline Distribution const & at(const Index &_index) const
    {

    }

    inline Distribution & at(const Index &_index)
    {

    }





private:
    Size             size;
    Resolution       resolution;
    Point            origin;
    std::size_t      data_size;
    Distribution    *data;



};
}
#endif // NDT_GRID_HPP
