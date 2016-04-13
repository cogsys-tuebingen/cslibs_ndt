#ifndef NDTGRID_H
#define NDTGRID_H

#include "laserscan.hpp"
#include "rolling_distribution.hpp"

namespace ndt {
class NDTGrid
{
public:
    typedef std::shared_ptr<NDTGrid> Ptr;

    NDTGrid(const std::size_t dim_x,
            const std::size_t dim_y,
            const double resolution,
            const double origin_off_x = 0.0,
            const double origin_off_y = 0.0);

    NDTGrid(const double size_x,
            const double size_y,
            const double resolution,
            const double origin_off_x = 0.0,
            const double origin_off_y = 0.0);

    inline void setOrigin(const Point &origin)
    {
        origin_[0] = origin(0);
        origin_[1] = origin(1);
    }

    inline bool add(const Point &point)
    {
        int idx = (point(0) - origin_[0]) / resolution_;
        int idy = (point(1) - origin_[1]) / resolution_;
        if(withinGrid(idx, idy)) {
            data_[idy * dim_x_ + idx].add(point);
            return true;
        }
        return false;
    }

    inline double sample(const Point &point)
    {
        std::size_t idx = (point(0) - origin_[0]) / resolution_;
        std::size_t idy = (point(1) - origin_[1]) / resolution_;
        if(withinGrid(idx, idy)) {
            return data_[idy * dim_x_ + idx].sample(point);
        }
        return 0.0;
    }

    inline bool sample(const Point &point,
                       double      &p,
                       Point       &mean,
                       Covariance  &inv_cov)
    {
        std::size_t idx = (point(0) - origin_[0]) / resolution_;
        std::size_t idy = (point(1) - origin_[1]) / resolution_;
        if(withinGrid(idx, idy)) {
            data_[idy * dim_x_ + idx].mean(mean);
            data_[idy * dim_x_ + idx].invCovariance(inv_cov);
            p = data_[idy * dim_x_ + idx].sample(point);
            return true;
        }
        return false;
    }

    inline void index(const Point &point,
                      std::size_t &idx,
                      std::size_t &idy)
    {
        idx = (point(0) - origin_[0]) / resolution_;
        idy = (point(1) - origin_[1]) / resolution_;
    }

    inline RollingDistribution & at (const std::size_t x,
                                     const std::size_t y)
    {
        return data_[y * dim_x_ + x];
    }

    inline RollingDistribution const & at (const std::size_t x,
                                           const std::size_t y) const
    {
        return data_[y * dim_x_ + x];
    }

    /**
     * @brief Check if point is within grid area or not.
     *        Attention: Be carefull with over or underflows!
     * @param idx - index x
     * @param idy - index y
     * @return  index in bounds
     */
    inline bool withinGrid(const int idx,
                           const int idy) const
    {
        return idx >= 0 &&
               idx < (int) dim_x_ &&
               idy >= 0 &&
               idy < (int) dim_y_;
    }

private:
    std::size_t dim_x_;
    std::size_t dim_y_;
    double      resolution_;
    double      origin_[2];

    RollingDistribution *data_;
};
}

#endif // NDTGRID_H
