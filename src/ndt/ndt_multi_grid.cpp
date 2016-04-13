#include "ndt_multi_grid.h"

using namespace ndt;
NDTMultiGrid::NDTMultiGrid(const std::size_t dim_x,
                           const std::size_t dim_y,
                           const double resolution) :
    data{NDTGrid(dim_x, dim_y, resolution),
         NDTGrid(dim_x, dim_y, resolution, 0.0, -resolution * 0.5),
         NDTGrid(dim_x, dim_y, resolution, -resolution * 0.5),
         NDTGrid(dim_x, dim_y, resolution, -resolution * 0.5, -resolution * 0.5)},
    resolution_(resolution)
{
}

NDTMultiGrid::NDTMultiGrid(const double size_x,
                           const double size_y,
                           const double resolution) :
    data{NDTGrid(size_x, size_y, resolution),
         NDTGrid(size_x, size_y, resolution, 0.0, -resolution * 0.5),
         NDTGrid(size_x, size_y, resolution, -resolution * 0.5),
         NDTGrid(size_x, size_y, resolution, -resolution * 0.5, -resolution * 0.5)},
    resolution_(resolution)
{
}
