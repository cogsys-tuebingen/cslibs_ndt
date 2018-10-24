#ifndef CSLIBS_NDT_3D_MATCH_DEFAULT_HPP
#define CSLIBS_NDT_3D_MATCH_DEFAULT_HPP

#include <cslibs_ndt_3d/matching/hessian.hpp>
#include <cslibs_ndt_3d/matching/jacobian.hpp>
#include <cslibs_ndt_3d/matching/params.hpp>
#include <cslibs_ndt_3d/matching/result.hpp>

#include <cslibs_math_3d/linear/pointcloud.hpp>
#include <cslibs_math_3d/algorithms/icp.hpp>

#include <cslibs_math/statistics/distribution.hpp>
#include <thread>
#include <atomic>

namespace cslibs_ndt_3d {
namespace matching {
namespace impl {

template<typename ndt_t, const std::size_t Ts>
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const typename ndt_t::ConstPtr               &dst,
                  const ParametersWithICP                      &params,
                  Result                                       &r)
{

}

template<typename ndt_t>
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const typename ndt_t::ConstPtr               &dst,
                  const ParametersWithICP                      &params,
                  Result                                       &r)
{

}


template<typename ndt_t, const std::size_t Ts>
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const typename ndt_t::ConstPtr               &dst,
                  const Parameters                             &params,
                  Result                                       &r)
{
    using point_t             = typename ndt_t::point_t;
    using point_transposed_t  = Eigen::Matrix<double, 1, 3>;
    using transform_t         = typename ndt_t::transform_t;
    using bundle_t            = typename ndt_t::distribution_bundle_t;
    using ndt_distribution_t  = typename ndt_t::distribution_t;
    using distribution_t      = typename ndt_distribution_t::distribution_t;
    using distributions_t     = std::array<distribution_t, 8>;
    using gradient_t          = Eigen::Matrix<double, 6, 1>;
    using hessian_t           = Eigen::Matrix<double, 6, 6>;

    auto get = [&dst](const point_t &p,
            distributions_t &d)
    {
        const bundle_t *b = dst->getDistributionBundle(p);
        if(b != nullptr) {
            for(std::size_t i = 0 ; i < 8 ; ++i) {
                d[i] = b->at(i)->getHandle()->data();
            }
            return true;
        }
        return false;
    };


    cslibs_math_3d::Pointcloud3d::Ptr src_prime(new cslibs_math_3d::Pointcloud3d(*src));
    src_prime->transform(params.transform());

    std::array<double,3> linear_old  = {{0.0, 0.0, 0.0}};
    std::array<double,3> angular_old = {{0.0, 0.0, 0.0}};

    std::array<double,3> linear_delta  = {{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}};
    std::array<double,3> angular_delta = {{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}};

    std::array<double,3> linear      = {{0.0, 0.0, 0.0}};
    std::array<double,3> angular     = {{0.0, 0.0, 0.0}};

    double max_score = std::numeric_limits<double>::lowest();
    double lambda = 1.0;

    std::size_t iterations = 0;
    std::size_t step_readjust = 0;
    Result::Termination termination = Result::ITERATIONS;

    /// termiantion tests
    const double rot_eps = params.rotEps();
    const double trans_eps = params.transEps();
    const std::size_t max_step_readjust = params.maxStepReadjust();
    const std::size_t max_iterations = params.maxIterations();
    const double alpha = params.alpha();

    auto eps = [&linear_delta, &angular_delta, rot_eps, trans_eps](){
        return linear_delta[0]  < trans_eps &&
                linear_delta[1]  < trans_eps &&
                linear_delta[2]  < trans_eps &&
                angular_delta[0] < rot_eps   &&
                angular_delta[1] < rot_eps   &&
                angular_delta[2] < rot_eps;

    };

    auto readjust = [&step_readjust, max_step_readjust]() {
        return step_readjust > 0 &&
                step_readjust > max_step_readjust;
    };

    const std::size_t size_prime = src_prime->size();

    for(std::size_t i = 0 ; i < max_iterations ; ++i, ++iterations) {
        if(readjust()) {
            /// drop out if too many consecutive step readjustments occur
            termination = Result::STEP_READJUSTMENTS;
            break;
        }

        Hessian  H;
        Hessian::get(angular, H);
        Jacobian J;
        Jacobian::get(angular, J);

        const transform_t t(linear[0],  linear[1],  linear[2],
                angular[0], angular[1], angular[2]);

        gradient_t g = gradient_t::Zero();
        hessian_t  h = hessian_t::Zero();
        gradient_t dp = gradient_t::Zero();
        double score = 0.0;

        std::atomic_uint  it_prime(0);

        std::mutex mutex_hg;
        auto exec = [size_prime, &it_prime, &mutex_hg, &h, &g, &score, &H, &J, &src_prime, &t, &get]()
        {
            std::size_t i = it_prime++;
            while(i < size_prime) {
                distributions_t distributions;
                const point_t x_prime = t * src_prime->at(i);
                gradient_t tg = gradient_t::Zero();
                hessian_t  th = hessian_t::Zero();
                double     ts = 0.0;

                if(get(x_prime, distributions)) {
                    for(const distribution_t &d : distributions) {
                        if(d.getN() < 3)
                            continue;

                        const typename distribution_t::sample_t            mu     = d.getMean();
                        const typename distribution_t::covariance_t        info   = d.getInformationMatrix();
                        const typename distribution_t::sample_t            q      = x_prime.data() - mu;
                        const typename distribution_t::sample_transposed_t q_info = static_cast<point_transposed_t>(q.transpose()) * info;
                        const double                                       e      = -0.5 * static_cast<double>(q_info * q);
                        const double                                       s      = std::exp(e);

                        for(std::size_t i = 0 ; i < 6 ; ++i) {
                            tg(i) += s * q_info * J.get(i, q);
                            for(std::size_t j = 0 ; j < 6 ; ++j) {
                                th(i,j) += s * static_cast<double>(q_info * H.get(i,j, q)) +
                                        s * static_cast<double>(static_cast<point_transposed_t>(J.get(i,q).transpose()) * info * J.get(j,q));
                            }
                        }
                        ts += s;
                    }
                    std::unique_lock<std::mutex> l(mutex_hg);
                    score += ts;
                    h += th;
                    g += tg;
                }
                i = it_prime++;
            }
        };

        std::array<std::thread, Ts> workers;
        for(std::thread &w : workers) {
            w = std::thread(exec);
        }
        for(std::thread &w : workers) {
            w.join();
        }

        if(score < max_score) {
            /// roll back and retry mit smaller step size
            lambda *= alpha;
            linear = linear_old;
            angular = angular_old;
            ++step_readjust;
            continue;
        }

        if(score > max_score) {
            /// update the max score and reduce step size
            max_score = score;
            lambda /= alpha;
            step_readjust = 0;
        }

        cslibs_math::statistics::LimitEigenValues<6, 3>::apply(h);
        dp  = -h.fullPivLu().solve(g);
        dp *= lambda;

        linear_old = linear;
        angular_old = angular;

        for(std::size_t i = 0 ; i < 3 ; ++i) {
            linear[i]  += dp(i);
            angular[i]  = cslibs_math::common::angle::normalize(angular[i] + dp(3 + i));
            linear_delta[i] = std::fabs(linear[i] - linear_old[i]);
            angular_delta[i] = std::fabs(cslibs_math::common::angle::difference(angular_old[i], angular[i]));
        }

        if(eps()) {
            termination = Result::EPS;
            break;
        }

    }

    const transform_t t(linear[0],  linear[1],  linear[2],
            angular[0], angular[1], angular[2]);

    r = Result(max_score, iterations, t * params.transform(), termination);
}

template<typename ndt_t, typename iterator_t>
inline void match(const iterator_t& src_begin,
                  const iterator_t& src_end,
                  const ndt_t&      dst,
                  const Parameters& params,
                  Result&           r)
{
    using point_t             = typename ndt_t::point_t;
    using point_transposed_t  = Eigen::Matrix<double, 1, 3>;
    using transform_t         = typename ndt_t::transform_t;
    using bundle_t            = typename ndt_t::distribution_bundle_t;
    using ndt_distribution_t  = typename ndt_t::distribution_t;
    using distribution_t      = typename ndt_distribution_t::distribution_t;
    using distributions_t     = std::array<distribution_t, 8>;
    using gradient_t          = Eigen::Matrix<double, 6, 1>;
    using hessian_t           = Eigen::Matrix<double, 6, 6>;

    auto get = [&dst](const point_t &p,
            distributions_t &d)
    {
        const bundle_t *b = dst.getDistributionBundle(p);
        if(b != nullptr) {
            for(std::size_t i = 0 ; i < 8 ; ++i) {
                d[i] = b->at(i)->getHandle()->data();
            }
            return true;
        }
        return false;
    };


    std::vector<point_t> src_prime;
    src_prime.reserve(std::distance(src_begin, src_end));
    std::transform(src_begin, src_end, std::back_inserter(src_prime),
                   [&](const point_t& point) { return params.transform() * point; });

    std::array<double,3> linear_old  = {{0.0, 0.0, 0.0}};
    std::array<double,3> angular_old = {{0.0, 0.0, 0.0}};

    std::array<double,3> linear_delta  = {{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}};
    std::array<double,3> angular_delta = {{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}};

    std::array<double,3> linear      = {{0.0, 0.0, 0.0}};
    std::array<double,3> angular     = {{0.0, 0.0, 0.0}};

    double max_score = std::numeric_limits<double>::lowest();
    double lambda = 1.0;

    std::size_t iterations = 0;
    std::size_t step_readjust = 0;
    Result::Termination termination = Result::ITERATIONS;

    /// termiantion tests
    const double rot_eps = params.rotEps();
    const double trans_eps = params.transEps();
    const std::size_t max_step_readjust = params.maxStepReadjust();
    const std::size_t max_iterations = params.maxIterations();
    const double alpha = params.alpha();

    auto eps = [&linear_delta, &angular_delta, rot_eps, trans_eps](){
        return linear_delta[0]  < trans_eps &&
                linear_delta[1]  < trans_eps &&
                linear_delta[2]  < trans_eps &&
                angular_delta[0] < rot_eps   &&
                angular_delta[1] < rot_eps   &&
                angular_delta[2] < rot_eps;

    };

    auto readjust = [&step_readjust, max_step_readjust]() {
        return step_readjust > 0 &&
                step_readjust > max_step_readjust;
    };

    for(std::size_t i = 0 ; i < max_iterations ; ++i, ++iterations) {
        if(readjust()) {
            /// drop out if too many consecutive step readjustments occur
            termination = Result::STEP_READJUSTMENTS;
            break;
        }

        Hessian  H;
        Hessian::get(angular, H);
        Jacobian J;
        Jacobian::get(angular, J);

        const transform_t t(linear[0],  linear[1],  linear[2],
                angular[0], angular[1], angular[2]);

        gradient_t g = gradient_t::Zero();
        hessian_t  h = hessian_t::Zero();

        gradient_t dp = gradient_t::Zero();

        distributions_t distributions;
        double score = 0.0;

        for(const point_t &x : src_prime) {
            const point_t x_prime = t * x;
            if(get(x_prime, distributions)) {
                for(const distribution_t &d : distributions) {
                    if(d.getN() < 3)
                        continue;

                    const typename distribution_t::sample_t            mu     = d.getMean();
                    const typename distribution_t::covariance_t        info   = d.getInformationMatrix();
                    const typename distribution_t::sample_t            q      = x_prime.data() - mu;
                    const typename distribution_t::sample_transposed_t q_info = static_cast<point_transposed_t>(q.transpose()) * info;
                    const double                                       e      = -0.5 * static_cast<double>(q_info * q);
                    const double                                       s      = std::exp(e);

                    for(std::size_t i = 0 ; i < 6 ; ++i) {
                        g(i) += s * q_info * J.get(i, q);
                        for(std::size_t j = 0 ; j < 6 ; ++j) {
                            h(i,j) += s * static_cast<double>(q_info * H.get(i,j, q)) +
                                    s * static_cast<double>(static_cast<point_transposed_t>(J.get(i,q).transpose()) * info * J.get(j,q));
                        }
                    }
                    score += s;
                }
            }
        }

        if(score < max_score) {
            /// roll back and retry mit smaller step size
            lambda *= alpha;
            linear = linear_old;
            angular = angular_old;
            ++step_readjust;
            continue;
        }

        if(score > max_score) {
            /// update the max score and reduce step size
            max_score = score;
            lambda /= alpha;
            step_readjust = 0;
        }

        cslibs_math::statistics::LimitEigenValues<6, 3>::apply(h);
        dp  = -h.fullPivLu().solve(g);
        dp *= lambda;

        linear_old = linear;
        angular_old = angular;

        for(std::size_t i = 0 ; i < 3 ; ++i) {
            linear[i]  += dp(i);
            angular[i]  = cslibs_math::common::angle::normalize(angular[i] + dp(3 + i));
            linear_delta[i] = std::fabs(linear[i] - linear_old[i]);
            angular_delta[i] = std::fabs(cslibs_math::common::angle::difference(angular_old[i], angular[i]));
        }

        if(eps()) {
            termination = Result::EPS;
            break;
        }

    }

    const transform_t t(linear[0],  linear[1],  linear[2],
            angular[0], angular[1], angular[2]);

    r = Result(max_score, iterations, t * params.transform(), termination);
}


template<typename ndt_t>
inline void match(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const typename ndt_t::ConstPtr               &dst,
                  const Parameters                             &params,
                  Result                                       &r)
{
    match(src->begin(), src->end(), *dst, params, r);
}

}
}
}

#endif // CSLIBS_NDT_3D_MATCH_DEFAULT_HPP
