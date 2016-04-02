#include "ndt_grid.h"

using namespace ndt;

NDTGrid::NDTGrid(const std::size_t dim_x,
                 const std::size_t dim_y,
                 const double resolution,
                 const double origin_off_x,
                 const double origin_off_y) :
    dim_x_(dim_x),
    dim_y_(dim_y),
    resolution_(resolution),
    origin_{-(dim_x_ * resolution * 0.5 + origin_off_x),
            -(dim_y_ * resolution * 0.5 + origin_off_y)},
  data_(new RollingDistribution[dim_x_ * dim_y_])
{
}

NDTGrid::NDTGrid(const double size_x,
                 const double size_y,
                 const double resolution,
                 const double origin_off_x,
                 const double origin_off_y) :
    dim_x_(size_x / resolution),
    dim_y_(size_y / resolution),
    resolution_(resolution),
    origin_{-(size_x * 0.5 + origin_off_x),
            -(size_y * 0.5 + origin_off_y)},
  data_(new RollingDistribution[dim_x_ * dim_y_])
{
}
