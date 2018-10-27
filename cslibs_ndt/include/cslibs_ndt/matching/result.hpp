#pragma once

#include <cstdint>
#include <eigen3/Eigen/Eigen>

namespace cslibs_ndt {
namespace matching {

enum class Termination { NONE, MAX_ITERATIONS, DELTA_EPSILON, MAX_STEP_READJUSTMENTS };

template<typename transform_t>
class EIGEN_ALIGN16 Result
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit Result() :
            Result(0.0, 0, transform_t{}, Termination::NONE)
    {}

    explicit Result(double score,
                    std::size_t iterations,
                    const transform_t& transform,
                    Termination termination) :
            score_(score),
            iterations_(iterations),
            transform_(transform),
            termination_(termination)
    {}

    double              score()         const { return score_; }
    std::size_t         iterations()    const { return iterations_; }
    const transform_t&  transform()     const { return transform_; }
    Termination         termination()   const { return termination_; }

    double&      score()        { return score_; }
    std::size_t& iterations()   { return iterations_; }
    transform_t& transform()    { return transform_; }
    Termination& termination()  { return termination_; }

protected:
    double      score_;
    std::size_t iterations_;
    transform_t transform_;
    Termination termination_;
};

}
}

namespace std {
std::string to_string(const cslibs_ndt::matching::Termination& termination)
{
    using namespace cslibs_ndt::matching;
    switch (termination)
    {
        default:
        case Termination::NONE: return "NONE";
        case Termination::MAX_ITERATIONS: return "MAX_ITERATIONS";
        case Termination::DELTA_EPSILON: return "DELTA_EPSILON";
        case Termination::MAX_STEP_READJUSTMENTS: return "MAX_STEP_READJUSTMENTS";
    }
}

template<typename transform_t>
std::string to_string(const cslibs_ndt::matching::Result<transform_t>& result)
{
    std::string s;
    s += "score      : " + std::to_string(result.score()) + "\n";
    s += "iterations : " + std::to_string(result.iterations()) + "\n";
    s += "transform  : " + std::to_string(result.transform()) + "\n";
    s += "termination: " + std::to_string(result.termination());
    return s;
}
}
