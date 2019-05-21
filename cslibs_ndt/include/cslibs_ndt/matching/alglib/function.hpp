#ifndef CSLIBS_NDT_MATCHING_ALGLIB_FUNCTION_HPP
#define CSLIBS_NDT_MATCHING_ALGLIB_FUNCTION_HPP

namespace cslibs_ndt {
namespace matching {
namespace alglib {

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
 *     };
 *
 *     // function calculating function value
 *     inline static void apply(const ::alglib::real_1d_array &x, double &f, void* ptr);
 *
 *     // function calculating function value and jacobian
 *     inline static void apply(const ::alglib::real_1d_array &x, double &f, ::alglib::real_1d_array &jac, void* ptr);
 *
 *     // function calculating function value, jacobian and hessian
 *     inline static void apply(const ::alglib::real_1d_array &x, double &f, ::alglib::real_1d_array &jac, ::alglib::real_2d_array &hes, void* ptr);
 * };
 */

}
}
}

#endif CSLIBS_NDT_MATCHING_ALGLIB_FUNCTION_HPP
