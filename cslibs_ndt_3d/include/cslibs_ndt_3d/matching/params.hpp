#ifndef CSLIBS_NDT_PARAMS_HPP
#define CSLIBS_NDT_PARAMS_HPP

#include <cslibs_math_3d/linear/transform.hpp>

namespace cslibs_ndt_3d {
namespace matching {
class EIGEN_ALIGN16 ParametersWithICP : public Parameters
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline ParametersWithICP(const double                        resolution = 0.5,
                             const std::size_t                   max_iterations = 100,
                             const std::size_t                   icp_max_iterations = 50,
                             const double                        icp_min_assigned_points = 0.75,
                             const double                        icp_max_distance = 0.5,
                             const double                        rot_eps = 1e-4,
                             const double                        trans_eps = 1e-4,
                             const std::size_t                   step_adjustment_retries = 5,
                             const double                        alpha = 1.0,
                             const cslibs_math_3d::Transform3d &transform = cslibs_math_3d::Transform3d()) :
        Parameters(resolution, max_iterations, rot_eps, trans_eps, step_adjustment_retries, alpha, transform),
        icp_max_iterations_(icp_max_iterations),
        icp_min_assigned_points_(icp_min_assigned_points),
        icp_max_distance_(icp_max_distance)
    {
    }

    inline std::size_t maxIterationsICP() const
    {
        return icp_max_iterations_;
    }

    inline std::size_t & maxIterationsICP()
    {
        return icp_max_iterations_;
    }

    inline double minAssignedPoints() const
    {
        return icp_min_assigned_points_;
    }

    inline double & minAssignedPoints()
    {
        return icp_min_assigned_points_;
    }

    inline double maxDistanceICP() const
    {
        return icp_max_distance_;
    }

    inline double & maxDistanceICP()
    {
        return icp_max_distance_;
    }

protected:
    std::size_t icp_max_iterations_;
    double      icp_min_assigned_points_;
    double      icp_max_distance_;
};
}
}

#endif // CSLIBS_NDT_PARAMS_HPP
