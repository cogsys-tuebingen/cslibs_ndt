#pragma once

namespace cslibs_ndt {
namespace matching {

template<typename MapT, typename Enable = void>
struct MatchTraits;
/*
Required Interface:
- "void"-usings have to be adjusted
- {LINEAR,ANGULAR}_DIMS have to be adjusted

{
    static constexpr int LINEAR_DIMS = 0;
    static constexpr int ANGULAR_DIMS = 0;
    using Jacobian = void;
    using Hessian  = void;

    using gradient_t = Eigen::Matrix<double, LINEAR_DIMS + ANGULAR_DIMS, 1>;
    using hessian_t  = Eigen::Matrix<double, LINEAR_DIMS + ANGULAR_DIMS, LINEAR_DIMS + ANGULAR_DIMS>;

    using point_t       = void;
    using transform_t   = void;

    static transform_t makeTransform(const Eigen::Matrix<double, LINEAR_DIMS, 1>& linear,
                                     const Eigen::Matrix<double, ANGULAR_DIMS, 1>& angular);

    static void computeGradient(const MapT& map,
                                const point_t& point,
                                const Jacobian& J,
                                const Hessian& H,
                                double& score,
                                gradient_t& g,
                                hessian_t& h);
};
*/
}
}
