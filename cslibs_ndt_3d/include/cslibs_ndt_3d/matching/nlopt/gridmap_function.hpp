#ifndef CSLIBS_NDT_3D_MATCHING_NLOPT_GRIDMAP_FUNCTION_HPP
#define CSLIBS_NDT_3D_MATCHING_NLOPT_GRIDMAP_FUNCTION_HPP

#include <cslibs_ndt/matching/nlopt/function.hpp>
#include <cslibs_ndt/map/map.hpp>

namespace cslibs_ndt {
namespace matching {
namespace nlopt {

template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t,
          typename point_t>
class Function<
        cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>,
        point_t> {
public:
    using ndt_t = cslibs_ndt::map::Map<option_t,3,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>;

    template <std::size_t n> // n=6: xyz rpy  (Euler);
    struct Functor {         // n=7: xyz xyzw (Quaternion)
        const ndt_t* map_;
        const std::vector<point_t>* points_;

        std::array<double,n> initial_guess_;
        double translation_weight_;
        double rotation_weight_;
        double map_weight_;
    };
    using FunctorRPY = Functor<6>;
    using FunctorQuaternion = Functor<7>;

    inline static double hypot(const double& x, const double& y, const double& z)
    {
        auto sq = [](const double& v) { return v*v; };
        return std::sqrt(sq(x) + sq(y) + sq(z));
    }

    inline static double applyRPY(unsigned n, const double *x, double *grad, void* ptr)
    {
        // since f is discontinuous, use derivative-free algorithms!
        if (grad) {
            std::cerr << "Gradient not implemented..." << std::endl;
            return 0;
        }

        // check that all necessary information is given
        const auto& casted_ptr = (Functor<6>*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return 0;
        }

        // dissolve Functor
        const Functor<6>& object = *casted_ptr;
        const auto& points       = *(object.points_);
        const auto& map          = *(object.map_);

        double fi = 0;
        const typename ndt_t::pose_t current_transform(x[0],x[1],x[2],x[3],x[4],x[5]); // xyz rpy

        // evaluate function
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1),p(2));
            const double score = map.sampleNonNormalized(q);
            fi += std::isnormal(score) ? (1.0 - score) : 1.0;
        }

        // calculate translational and rotational function component
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1], x[2] - initial_guess[2]); // xyz
        const double rot_diff     = hypot(cslibs_math::common::angle::difference(x[3], initial_guess[3]),  // rpy
                                          cslibs_math::common::angle::difference(x[4], initial_guess[4]),
                                          cslibs_math::common::angle::difference(x[5], initial_guess[5]));

        // apply weights
        fi = object.map_weight_ * fi / static_cast<double>(points.size()) +
             object.translation_weight_ * trans_diff +
             object.rotation_weight_ * std::fabs(rot_diff);

        return fi;
    }

    inline static double mapScoreRPY(const double *x, const double &fi, void* ptr)
    {
        // check that all necessary information is given
        const auto& casted_ptr = (Functor<6>*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return 0;
        }

        // calculate translational and rotational function component
        const Functor<6>& object  = *casted_ptr;
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1], x[2] - initial_guess[2]); // xyz
        const double rot_diff     = hypot(cslibs_math::common::angle::difference(x[3], initial_guess[3]),  // rpy
                                          cslibs_math::common::angle::difference(x[4], initial_guess[4]),
                                          cslibs_math::common::angle::difference(x[5], initial_guess[5]));

        // extract real fvalue (map correlation value)
        return 1.0 - (fi -
                object.translation_weight_ * trans_diff -
                object.rotation_weight_ * std::fabs(rot_diff)) /
                object.map_weight_;
    }

    inline static double applyQuaternion(unsigned n, const double *x, double *grad, void* ptr)
    {
        // since f is discontinuous, use derivative-free algorithms!
        if (grad) {
            std::cerr << "Gradient not implemented..." << std::endl;
            return 0;
        }

        // check that all necessary information is given
        const auto& casted_ptr = (Functor<7>*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return 0;
        }

        // dissolve Functor
        const Functor<7>& object = *casted_ptr;
        const auto& points       = *(object.points_);
        const auto& map          = *(object.map_);

        double fi = 0;
        const typename ndt_t::pose_t current_transform(
                    cslibs_math_3d::Vector3<_T>(x[0],x[1],x[2]),          // xyz
                    cslibs_math_3d::Quaternion<_T>(x[3],x[4],x[5],x[6])); // xyzw
        const auto& rot = current_transform.rotation();

        // evaluate function
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1),p(2));
            const double score = map.sampleNonNormalized(q);
            fi += std::isnormal(score) ? (1.0 - score) : 1.0;
        }

        // calculate translational and rotational function component
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1], x[2] - initial_guess[2]); // xyz
        const double rot_diff     = hypot(cslibs_math::common::angle::difference(rot.roll(), initial_guess[3]),  // rpy
                                          cslibs_math::common::angle::difference(rot.pitch(), initial_guess[4]),
                                          cslibs_math::common::angle::difference(rot.yaw(), initial_guess[5]));

        // apply weights
        fi = object.map_weight_ * fi / static_cast<double>(points.size()) +
             object.translation_weight_ * trans_diff +
             object.rotation_weight_ * std::fabs(rot_diff);

        return fi;
    }

    inline static double mapScoreQuaternion(const double *x, const double &fi, void* ptr)
    {
        // check that all necessary information is given
        const auto& casted_ptr = (Functor<7>*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return 0;
        }
        const cslibs_math_3d::Quaternion<_T> rot(x[3],x[4],x[5],x[6]);

        // calculate translational and rotational function component
        const Functor<7>& object  = *casted_ptr;
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1], x[2] - initial_guess[2]); // xyz
        const double rot_diff     = hypot(cslibs_math::common::angle::difference(rot.roll(), initial_guess[3]),  // rpy
                                          cslibs_math::common::angle::difference(rot.pitch(), initial_guess[4]),
                                          cslibs_math::common::angle::difference(rot.yaw(), initial_guess[5]));

        // extract real fvalue (map correlation value)
        return 1.0 - (fi -
                object.translation_weight_ * trans_diff -
                object.rotation_weight_ * std::fabs(rot_diff)) /
                object.map_weight_;
    }
};

}
}
}

#endif // CSLIBS_NDT_3D_MATCHING_NLOPT_GRIDMAP_FUNCTION_HPP
