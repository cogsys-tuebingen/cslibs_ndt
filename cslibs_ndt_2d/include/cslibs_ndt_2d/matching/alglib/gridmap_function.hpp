#ifndef CSLIBS_NDT_2D_MATCHING_ALGLIB_GRIDMAP_FUNCTION_HPP
#define CSLIBS_NDT_2D_MATCHING_ALGLIB_GRIDMAP_FUNCTION_HPP

#include <cslibs_ndt/matching/alglib/function.hpp>
#include <cslibs_ndt/map/map.hpp>

#include <optimization.h>

namespace cslibs_ndt {
namespace matching {
namespace alglib {

template <cslibs_ndt::map::tags::option option_t,
          typename _T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t,
          typename point_t>
class Function<
        cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>,
        point_t> {
public:
    using ndt_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::Distribution,_T,backend_t,dynamic_backend_t>;

    struct Functor {
        const ndt_t* map_;
        const std::vector<point_t>* points_;

        std::array<double,3> initial_guess_;
        double translation_weight_;
        double rotation_weight_;
        double map_weight_;
    };

    inline static void apply(const ::alglib::real_1d_array &x, ::alglib::real_1d_array &fi, void *ptr)
    {
        // check that all necessary information is given
        const auto& casted_ptr = (Functor*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return;
        }

        // dissolve Functor
        const Functor& object = *casted_ptr;
        const auto& points    = *(object.points_);
        const auto& map       = *(object.map_);

        fi[0] = 0;
        const typename ndt_t::pose_t current_transform(x[0],x[1],x[2]);

        // evaluate function
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1));
            const double score = map.sampleNonNormalized(q);
            fi[0] += std::isnormal(score) ? (1.0 - score) : 1.0;
        }

        // calculate translational and rotational function component
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1]);
        const double rot_diff     = cslibs_math::common::angle::difference(x[2], initial_guess[2]);

        // apply weights
        fi[0] = object.map_weight_ * fi[0] / static_cast<double>(points.size()) +
                object.translation_weight_ * trans_diff +
                object.rotation_weight_ * std::fabs(rot_diff);
        fi[0] = std::sqrt(std::fabs(fi[0]));
    }

    inline static double mapScore(const ::alglib::real_1d_array &x, const ::alglib::real_1d_array &fi, void* ptr)
    {
        // check that all necessary information is given
        const auto& casted_ptr = (Functor*)ptr;
        if (!casted_ptr) {
            std::cerr << "Correct Functor not given..." << std::endl;
            return 0;
        }

        // calculate translational and rotational function component
        const Functor& object     = *casted_ptr;
        const auto& initial_guess = (object.initial_guess_);
        const double trans_diff   = hypot(x[0] - initial_guess[0], x[1] - initial_guess[1]);
        const double rot_diff     = cslibs_math::common::angle::difference(x[2], initial_guess[2]);

        // extract real fvalue (map correlation value)
        return 1.0 - (fi[0] -
                object.translation_weight_ * trans_diff -
                object.rotation_weight_ * std::fabs(rot_diff)) /
                object.map_weight_;
    }
};

}
}
}

#endif // CSLIBS_NDT_2D_MATCHING_ALGLIB_GRIDMAP_FUNCTION_HPP
