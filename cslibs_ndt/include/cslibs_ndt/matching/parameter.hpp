#pragma once

namespace cslibs_ndt {
namespace matching {

struct Parameter
{
    std::size_t max_iterations;
    double translation_epsilon;
    double rotation_epsilon;
    std::size_t max_step_readjustments;
    double alpha;
};

}
}
