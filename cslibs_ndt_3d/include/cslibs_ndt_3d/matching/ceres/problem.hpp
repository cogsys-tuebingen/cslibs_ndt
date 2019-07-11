#include <cslibs_ndt/matching/ceres/problem.hpp>

#include <cslibs_ndt_3d/matching/ceres/map/gridmap_cost_functor.hpp>
#include <cslibs_ndt_3d/matching/ceres/map/occupancy_gridmap_cost_functor.hpp>

#include <ceres/problem.h>
#include <ceres/solver.h>
#include <ceres/covariance.h>
