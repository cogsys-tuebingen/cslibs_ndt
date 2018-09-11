#ifndef CSLIBS_NDT_RESULT_HPP
#define CSLIBS_NDT_RESULT_HPP

#include <cslibs_math_3d/linear/transform.hpp>

namespace cslibs_ndt_3d {
namespace matching {
class EIGEN_ALIGN16 Result {
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    enum Termination {eps, iteration, step_readjust};

    inline Result() :
        score_(0.0),
        iterations_(0),
        termination_(iteration)
    {
    }

    inline Result(const double score,
                  const std::size_t iterations,
                  const cslibs_math_3d::Transform3d &transform,
                  const Termination termination) :
        score_(score),
        iterations_(iterations),
        transform_(transform),
        termination_(termination)
    {
    }

    inline double getScore() const
    {
        return score_;
    }

    inline std::size_t getIterations() const
    {
        return iterations_;
    }

    inline const cslibs_math_3d::Transform3d& getTransform() const
    {
        return transform_;
    }

    inline Termination getTermination() const
    {
        return termination_;
    }

    inline void setScore(const double s)
    {
        score_ = s;
    }

    inline void setIterations(const std::size_t i)
    {
        iterations_ = i;
    }

    inline void setTransform(const cslibs_math_3d::Transform3d &t)
    {
        transform_ = t;
    }

    inline void setTermination(const Termination &t)
    {
        termination_ = t;
    }

private:
    double                      score_;
    std::size_t                 iterations_;
    cslibs_math_3d::Transform3d transform_;
    Termination                 termination_;

};
}
}

namespace std {
inline std::string to_string(const cslibs_ndt_3d::matching::Result &r)
{
    std::string s;
    s += "score       " + std::to_string(r.getScore())      + "\n";
    s += "iterations  " + std::to_string(r.getIterations()) + "\n";
    s += "transform   " + std::to_string(r.getTransform())  + "\n";
    s += "termiantion ";
    switch(r.getTermination()) {
    case cslibs_ndt_3d::matching::Result::eps:
        s += "eps";
        break;
    case cslibs_ndt_3d::matching::Result::iteration:
        s += "iteration";
        break;
    case cslibs_ndt_3d::matching::Result::step_readjust:
        s+= "step_readjust";
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
}

#endif // CSLIBS_NDT_RESULT_HPP
