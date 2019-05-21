#ifndef CSLIBS_NDT_2D_MATCHING_ALGLIB_OCCUPANCY_GRIDMAP_FUNCTION_HPP
#define CSLIBS_NDT_2D_MATCHING_ALGLIB_OCCUPANCY_GRIDMAP_FUNCTION_HPP

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
        cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>,
        point_t> {
public:
    using ndt_t = cslibs_ndt::map::Map<option_t,2,cslibs_ndt::OccupancyDistribution,_T,backend_t,dynamic_backend_t>;

    struct Functor {
        const ndt_t* map_;
        const std::vector<point_t>* points_;
        const typename ndt_t::inverse_sensor_model_t::Ptr* ivm_;
    };

    inline static void apply(const ::alglib::real_1d_array &x, double &f, void* ptr)
    {
        const Functor& object = *((Functor*)ptr);
        const auto& points    = *(object.points_);
        const auto& map       = *(object.map_);
        const auto& ivm       = *(object.ivm_);
        const auto& orig_inv  = map.getInitialOrigin().inverse();

        f = 0;

        const typename ndt_t::pose_t current_transform(x[0],x[1],x[2]);
        std::size_t count = 0;
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1));

            const auto& bundle = map.get(q);
            if (!bundle)
                continue;
            ++count;

            const typename ndt_t::point_t q_prime = orig_inv * q;
            for (std::size_t i=0; i<4; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const auto& mean = di->getMean();
                        const auto inf = di->getInformationMatrix();

                        const double occ = bi->getOccupancy(ivm);
                        const auto& qm = q_prime.data() - mean;
                        const auto& qm_inf = qm.transpose() * inf;
                        const auto& exponent = -0.5 * qm_inf * qm;
                        const double score = 0.25 * std::exp(exponent) * occ;
                        f -= score;
                    }
                }
            }
        }

        f /= static_cast<double>(count);
    }

    inline static void apply(const ::alglib::real_1d_array &x, double &f, ::alglib::real_1d_array &jac, void* ptr)
    {
        const Functor& object = *((Functor*)ptr);
        const auto& points    = *(object.points_);
        const auto& map       = *(object.map_);
        const auto& ivm       = *(object.ivm_);
        const auto& orig_inv  = map.getInitialOrigin().inverse();

        f = 0;
        for (int i=0; i<3; ++i)
          jac[i] = 0;

        const typename ndt_t::pose_t current_transform(x[0],x[1],x[2]);
        const auto& s = current_transform.sin();
        const auto& c = current_transform.cos();

        std::size_t count = 0;
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1));

            const auto& bundle = map.get(q);
            if (!bundle)
                continue;
            ++count;

            const typename ndt_t::point_t q_prime = orig_inv * q;
            for (std::size_t i=0; i<4; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const auto& mean = di->getMean();
                        const auto inf = di->getInformationMatrix();

                        const double occ = bi->getOccupancy(ivm);
                        const auto& qm = q_prime.data() - mean;
                        const auto& qm_inf = qm.transpose() * inf;
                        const auto& exponent = -0.5 * qm_inf * qm;
                        const double score = 0.25 * std::exp(exponent) * occ;
                        f -= score;

                        std::array<Eigen::Matrix<_T,2,1>,3> d;
                        d[0] << 1, 0;
                        d[1] << 0, 1;
                        d[2] << -s*qm(0)-c*qm(1), c*qm(0)-s*qm(1);
                        jac[0] += score * qm_inf * d[0];
                        jac[1] += score * qm_inf * d[1];
                        jac[2] += score * qm_inf * d[2];
                    }
                }
            }
        }

        f /= static_cast<double>(count);
        for (int i=0; i<3; ++i)
          jac[i] /= static_cast<double>(count);
    }

    inline static void apply(const ::alglib::real_1d_array &x, double &f, ::alglib::real_1d_array &jac, ::alglib::real_2d_array &hes, void* ptr)
    {
        const Functor& object = *((Functor*)ptr);
        const auto& points    = *(object.points_);
        const auto& map       = *(object.map_);
        const auto& ivm       = *(object.ivm_);
        const auto& orig_inv  = map.getInitialOrigin().inverse();

        f = 0;
        for (int i=0; i<3; ++i) {
          jac[i] = 0;
          for (int j=0; j<3; ++j)
            hes[i][j] = 0;
        }

        const typename ndt_t::pose_t current_transform(x[0],x[1],x[2]);
        const auto& s = current_transform.sin();
        const auto& c = current_transform.cos();

        std::size_t count = 0;
        for (const auto& p : points) {
            const typename ndt_t::point_t q = current_transform * typename ndt_t::point_t(p(0),p(1));

            const auto& bundle = map.get(q);
            if (!bundle)
                continue;
            ++count;

            const typename ndt_t::point_t q_prime = orig_inv * q;
            for (std::size_t i=0; i<4; ++i) {
                if (const auto& bi = bundle->at(i)) {
                    if (const auto& di = bi->getDistribution()) {
                        if (!di->valid())
                            continue;

                        const auto& mean = di->getMean();
                        const auto inf = di->getInformationMatrix();

                        const double occ = bi->getOccupancy(ivm);
                        const auto& qm = q_prime.data() - mean;
                        const auto& qm_inf = qm.transpose() * inf;
                        const auto& exponent = -0.5 * qm_inf * qm;
                        const double score = 0.25 * std::exp(exponent) * occ;
                        f -= score;

                        std::array<Eigen::Matrix<_T,2,1>,3> d;
                        d[0] << 1, 0;
                        d[1] << 0, 1;
                        d[2] << -s*qm(0)-c*qm(1), c*qm(0)-s*qm(1);
                        jac[0] += score * qm_inf * d[0];
                        jac[1] += score * qm_inf * d[1];
                        jac[2] += score * qm_inf * d[2];

                        Eigen::Matrix<_T,2,1> h;
                        h << -c*qm[0]+s*qm[1], -s*qm[0]-c*qm[1];
                        hes[0][0] += score * ((qm_inf*d[0]).value()*(-qm_inf*d[0]).value() + d[0].transpose()*inf*d[0]);
                        hes[0][1] += score * ((qm_inf*d[0]).value()*(-qm_inf*d[1]).value() + d[1].transpose()*inf*d[0]);
                        hes[0][2] += score * ((qm_inf*d[0]).value()*(-qm_inf*d[2]).value() + d[2].transpose()*inf*d[0]);
                        hes[1][0] += score * ((qm_inf*d[1]).value()*(-qm_inf*d[0]).value() + d[0].transpose()*inf*d[1]);
                        hes[1][1] += score * ((qm_inf*d[1]).value()*(-qm_inf*d[1]).value() + d[1].transpose()*inf*d[1]);
                        hes[1][2] += score * ((qm_inf*d[1]).value()*(-qm_inf*d[2]).value() + d[2].transpose()*inf*d[1]);
                        hes[2][0] += score * ((qm_inf*d[2]).value()*(-qm_inf*d[0]).value() + d[0].transpose()*inf*d[2]);
                        hes[2][1] += score * ((qm_inf*d[2]).value()*(-qm_inf*d[1]).value() + d[1].transpose()*inf*d[2]);
                        hes[2][2] += score * ((qm_inf*d[2]).value()*(-qm_inf*d[2]).value() + qm_inf*h + d[2].transpose()*inf*d[2]);
                    }
                }
            }
        }

        f /= static_cast<double>(count);
        for (int i=0; i<3; ++i) {
          jac[i] /= static_cast<double>(count);
          for (int j=0; j<3; ++j)
            hes[i][j] /= static_cast<double>(count);
        }
    }
};

}
}
}

#endif // CSLIBS_NDT_2D_MATCHING_ALGLIB_OCCUPANCY_GRIDMAP_FUNCTION_HPP
