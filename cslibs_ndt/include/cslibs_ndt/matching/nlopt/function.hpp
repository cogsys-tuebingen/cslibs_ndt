#ifndef CSLIBS_NDT_MATCHING_NLOPT_FUNCTION_HPP
#define CSLIBS_NDT_MATCHING_NLOPT_FUNCTION_HPP

namespace cslibs_ndt {
namespace matching {
namespace nlopt {

template <typename ndt_t, typename point_t>
class Function;
/*
 * // Required Interface:
 * {
 * public:
 *     // helper object for all necessary information
 *     struct Functor {
 *         const ndt_t* map_;
 *         const std::vector<point_t>* points_;
 *         ...
 *     };
 *
 *     // function calculating function value (and gradient, if grad != NULL)
 *     // x:      current transformation estimate
 *     // grad:   container for multi-dimensional gradient
 *     // return: function value
 *     inline static double apply(unsigned n, const double *x, double *grad, void* ptr);
 * };
 */

}
}
}

// partial specializations can be found in
// cslibs_ndt_2d,
// cslibs_ndt_3d

#endif // CSLIBS_NDT_MATCHING_NLOPT_FUNCTION_HPP
