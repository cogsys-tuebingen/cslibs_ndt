#ifndef CSLIBS_NDT_RESULT_HPP
#define CSLIBS_NDT_RESULT_HPP

#include <cslibs_math_3d/linear/transform.hpp>

#include <cslibs_ndt/matching/result.hpp>

namespace cslibs_ndt_3d {
namespace matching {

enum class ICPTermination {NONE, MAX_ITERATIONS, DELTA_EPS, ASSIGNMENT_SUCCESS};

class EIGEN_ALIGN16 ResultWithICP : public cslibs_ndt::matching::Result<cslibs_math_3d::Transform3d>
{
public:
    using base_t = typename cslibs_ndt::matching::Result<cslibs_math_3d::Transform3d>;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline explicit ResultWithICP() :
        Result(),
        icp_iterations_(0),
        icp_termination_(ICPTermination::NONE)
    {
    }

    inline explicit ResultWithICP(const double                              score,
                                  const std::size_t                         iterations,
                                  const cslibs_math_3d::Transform3d        &transform,
                                  const cslibs_ndt::matching::Termination   termination,
                                  const std::size_t                         icp_iterations,
                                  const cslibs_math_3d::Transform3d        &icp_transform,
                                  const ICPTermination                      icp_termination) :
        base_t(score, iterations, transform, termination),
        icp_transform_(icp_transform),
        icp_iterations_(icp_iterations),
        icp_termination_(icp_termination)
    {
    }

    inline void assign(const base_t &b)
    {
        score_          = b.score();
        iterations_     = b.iterations();
        transform_      = b.transform();
        termination_    = b.termination();
    }

    inline std::size_t & icpIterations()
    {
        return icp_iterations_;
    }

    inline std::size_t icpIterations() const
    {
        return icp_iterations_;
    }

    inline ICPTermination & icpTermination()
    {
        return icp_termination_;
    }

    inline ICPTermination icpTermination() const
    {
        return icp_termination_;
    }

    inline cslibs_math_3d::Transform3d & ICPTransform()
    {
        return icp_transform_;
    }

    inline const cslibs_math_3d::Transform3d& ICPTransform() const
    {
        return icp_transform_;
    }

    inline Eigen::Matrix3d& icpCovariance()
    {
        return icp_covariance_;
    }

    inline const Eigen::Matrix3d& icpCovariance() const
    {
        return icp_covariance_;
    }

protected:
    cslibs_math_3d::Transform3d     icp_transform_;
    std::size_t                     icp_iterations_;
    ICPTermination                  icp_termination_;
    Eigen::Matrix3d                 icp_covariance_;

};
}
}

#endif // CSLIBS_NDT_RESULT_HPP
