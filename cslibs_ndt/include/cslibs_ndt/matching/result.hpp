#pragma once

namespace cslibs_ndt {
namespace matching {

enum class Termination { MAX_ITERATIONS, DELTA_EPSILON, MAX_STEP_READJUSTMENTS };

template<typename transform_t>
struct Result
{
    double score;
    std::size_t iterations;
    transform_t transform;
    Termination termination;
};

}
}
