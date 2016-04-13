#ifndef NDTMULTIGRID_H
#define NDTMULTIGRID_H

#include "ndt_grid.h"

namespace ndt {
class NDTMultiGrid
{
public:
    typedef std::shared_ptr<NDTMultiGrid> Ptr;

    NDTMultiGrid(const std::size_t dim_x,
                 const std::size_t dim_y,
                 const double resolution);

    NDTMultiGrid(const double size_x,
                 const double size_y,
                 const double resolution);

    inline bool add(const Point &p)
    {
        return  data[0].add(p) ||
                data[1].add(p) ||
                data[2].add(p) ||
                data[3].add(p);
    }

    inline double sample(const Point &point)
    {
        return data[0].sample(point) +
               data[1].sample(point) +
               data[2].sample(point) +
               data[3].sample(point);
    }

    NDTGrid & at(const std::size_t i)
    {
        return data[i];
    }

    NDTGrid & leftUpper()
    {
        return data[0];
    }

    NDTGrid const & leftUpper() const
    {
        return data[0];
    }

    NDTGrid & leftLower()
    {
        return data[1];
    }

    NDTGrid const & leftLower() const
    {
        return data[1];
    }

    NDTGrid & rightUpper()
    {
        return data[2];
    }

    NDTGrid const & rightUpper() const
    {
        return data[2];
    }

    NDTGrid & rightLower()
    {
        return data[3];
    }

    NDTGrid const & rightLower() const
    {
        return data[3];
    }

    inline double resolution() const
    {
        return resolution_;
    }

private:
    NDTGrid data[4];

    double resolution_;


};
}

#endif // NDTMULTIGRID_H
