#pragma once

#include <cstdint>

namespace cslibs_ndt {
namespace matching {

class Parameter
{
public:
    Parameter() :
        max_iterations_(35),
        translation_epsilon_(1e-3),
        rotation_epsilon_(1e-3),
        max_step_readjustments_(5),
        alpha_(1.1)
    {
    }

    explicit Parameter(std::size_t max_iterations,
                       double translation_epsilon,
                       double rotation_epsilon,
                       std::size_t max_step_readjustments,
                       double alpha) :
            max_iterations_(max_iterations),
            translation_epsilon_(translation_epsilon),
            rotation_epsilon_(rotation_epsilon),
            max_step_readjustments_(max_step_readjustments),
            alpha_(alpha)
    {}

    std::size_t maxIterations() const { return max_iterations_; }
    double translationEpsilon() const { return translation_epsilon_; }
    double rotationEpsilon() const { return rotation_epsilon_; }
    std::size_t maxStepReadjustments() const { return max_step_readjustments_; }
    double alpha() const { return alpha_; }

    std::size_t& maxIterations() { return max_iterations_; }
    double& translationEpsilon() { return translation_epsilon_; }
    double& rotationEpsilon() { return rotation_epsilon_; }
    std::size_t& maxStepReadjustments() { return max_step_readjustments_; }
    double& alpha() { return alpha_; }


private:
    std::size_t max_iterations_;
    double translation_epsilon_;
    double rotation_epsilon_;
    std::size_t max_step_readjustments_;
    double alpha_;
};

}
}
