#ifndef CSLIBS_NDT_RESULT_HPP
#define CSLIBS_NDT_RESULT_HPP

#include <cslibs_math_3d/linear/transform.hpp>

namespace cslibs_ndt_3d {
namespace matching {
class EIGEN_ALIGN16 ResultWithICP : public Result {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    enum ICPTermination {ICP_EPS, ICP_ITERATIONS, ICP_ASSIGNMENT};

    inline ResultWithICP() :
        Result(),
        icp_iterations_(0),
        icp_termination_(ICP_ITERATIONS)
    {
    }

    inline ResultWithICP(const double                       score,
                         const std::size_t                  iterations,
                         const cslibs_math_3d::Transform3d &transform,
                         const Result::Termination          termination,
                         const std::size_t                  icp_iterations,
                         const cslibs_math_3d::Transform3d &icp_transform,
                         const ICPTermination               icp_termination) :
        Result(score, iterations, transform, termination),
        icp_transform_(icp_transform),
        icp_iterations_(icp_iterations),
        icp_termination_(icp_termination)
    {
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

namespace std {
inline std::string to_string(const cslibs_ndt_3d::matching::Result &r)
{
    std::string s;
    s += "score          " + std::to_string(r.score())      + "\n";
    s += "iterations     " + std::to_string(r.iterations()) + "\n";
    s += "transform      " + std::to_string(r.transform())  + "\n";
    s += "termination    ";
    switch(r.termination()) {
    case cslibs_ndt_3d::matching::Result::EPS:
        s += "EPS";
        break;
    case cslibs_ndt_3d::matching::Result::ITERATIONS:
        s += "ITERATIONS";
        break;
    case cslibs_ndt_3d::matching::Result::STEP_READJUSTMENTS:
        s+= "STEP_READJUSTMENTS";
        break;
    default:
        break;
    }
    s += "\n";
    return s;
}

inline std::string to_string(const cslibs_ndt_3d::matching::ResultWithICP &r)
{
    std::string s;
    s += to_string(static_cast<const cslibs_ndt_3d::matching::Result&>(r));
    s += "icp iterations  " + std::to_string(r.icpIterations()) + "\n";
    s += "icp termination ";
    switch(r.icpTermination()) {
    case cslibs_ndt_3d::matching::ResultWithICP::ICP_ITERATIONS:
        s += "ITERATIONS";
        break;
    case cslibs_ndt_3d::matching::ResultWithICP::ICP_EPS:
        s += "EPS";
        break;
    case cslibs_ndt_3d::matching::ResultWithICP::ICP_ASSIGNMENT:
        s+= "ERROR";
        break;
    default:
        break;
    }
    s += "\n";
    return s;
}

inline std::ostream & operator << (std::ostream &out, const cslibs_ndt_3d::matching::Result &r)
{
    out << std::to_string(r);
    return out;
}

inline std::ostream & operator << (std::ostream &out, const cslibs_ndt_3d::matching::ResultWithICP &r)
{
    out << std::to_string(r);
    return out;
}
}

#endif // CSLIBS_NDT_RESULT_HPP
