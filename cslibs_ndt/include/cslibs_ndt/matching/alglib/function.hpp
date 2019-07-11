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
 *         ...
 *     };
 *
 *     // function calculating function value
 *     inline static void apply(const ::alglib::real_1d_array &x, ::alglib::real_1d_array &fi, void *ptr);
 * };
 */

}
}
}

// partial specializations can be found in
// cslibs_ndt_2d,
// cslibs_ndt_3d

#endif // CSLIBS_NDT_MATCHING_ALGLIB_FUNCTION_HPP
