#ifndef MATCHER_HPP
#define MATCHER_HPP
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>

#include <ndt/data/pointcloud.hpp>

namespace ndt {
namespace matching {
template<std::size_t Dim>
class Matcher {
public:
    typedef std::shared_ptr<Matcher>                     Ptr;
    typedef data::Pointcloud<Dim>                        PointCloudType;
    typedef std::array<double, Dim>                      ResolutionType;
    typedef Eigen::Matrix<double, 6, 1>                  LambdaType;        /// TODO : Move lambda to matcher directlyht
    typedef Eigen::Transform<double, Dim, Eigen::Affine> TransformType;
    typedef Matcher<Dim>                                 BaseClass;

    struct Parameters {
        Parameters() :
            eps_rot(1e-3),
            eps_trans(1e-3),
            max_iterations(100),
            max_step_corrections(10),
            lambda(LambdaType::Constant(0.1)),
            alpha(2.0)
        {
            resolution.fill(1.0);
        }
        Parameters(const ResolutionType &_resolution,
                   const double          _eps_rot,
                   const double          _eps_trans,
                   const std::size_t     _max_iterations,
                   const std::size_t     _max_step_corrections,
                   const LambdaType &_lambda,
                   const double _alpha) :
            resolution(_resolution),
            eps_rot(_eps_rot),
            eps_trans(_eps_trans),
            max_iterations(_max_iterations),
            max_step_corrections(_max_step_corrections),
            lambda(_lambda),
            alpha(_alpha)
        {
        }

        ResolutionType          resolution;
        double                  eps_rot;
        double                  eps_trans;
        std::size_t             max_iterations;
        std::size_t             max_step_corrections;
        LambdaType              lambda;
        double                  alpha;

    };



    Matcher(const Parameters &_params = Parameters()) :
        params(_params)
    {
    }

    /**
     * @brief match
     * @param _dst - point set to match onto
     * @param _src - point set to be matched onto _dst
     * @param _transformation - transformation from src frame to dst frame
     * @param _prior_transformation - initial transformation
     * @return
     */
    inline virtual double match(const PointCloudType &_dst,
                                const PointCloudType &_src,
                                TransformType        &_transformation,
                                const TransformType  &_prior_transformation = TransformType::Identity()) = 0;

protected:
    Parameters                params;

    inline bool epsTrans(const double a,
                         const double b) const
    {
        return fabs(a - b) < params.eps_trans;
    }

    inline bool epsRot(const double a,
                       const double b) const
    {
        return fabs(a - b) < params.eps_rot;
    }

    inline bool eps(const double a,
                    const double b,
                    const double e)
    {
        return fabs(a - b) < e;
    }

};
}
}
#endif // MATCHER_HPP
