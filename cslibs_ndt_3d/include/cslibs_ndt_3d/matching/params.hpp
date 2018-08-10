#ifndef CSLIBS_NDT_PARAMS_HPP
#define CSLIBS_NDT_PARAMS_HPP

#include <cslibs_math_3d/linear/transform.hpp>

namespace cslibs_ndt_3d {
namespace matching {
class Parameters {
public:
  inline Parameters(const double resolution_ = 0.5,
                   const std::size_t max_iterations_ = 100,
                   const double rot_eps_ = 1e-4,
                   const double trans_eps_ = 1e-4,
                   const std::size_t step_adjustment_retries = 5,
                   const double alpha_ = 1.0,
                   const cslibs_math_3d::Transform3d &transform_ = cslibs_math_3d::Transform3d()) :
    resolution_(resolution_),
    max_iterations_(max_iterations_),
    rot_eps_(rot_eps_),
    trans_eps_(trans_eps_),
    max_step_readjust_(step_adjustment_retries),
    alpha_(alpha_),
    transform_(transform_)
  {
  }

  inline double getResolution() const
  {
    return resolution_;
  }

  inline std::size_t getMaxIterations() const
  {
    return max_iterations_;
  }

  inline double getRotEps() const
  {
    return rot_eps_;
  }

  inline double getTransEps() const
  {
    return trans_eps_;
  }

  inline std::size_t getMaxStepReadjust() const
  {
    return max_step_readjust_;
  }

  inline double getAlpha() const
  {
      return alpha_;
  }

  inline const cslibs_math_3d::Transform3d &getTransform() const
  {
      return transform_;
  }

  inline void setResolution(const double r)
  {
    resolution_ = r;
  }

  inline void setMaxIterations(const std::size_t i)
  {
    max_iterations_ = i;
  }

  inline void setRotEps(const double r)
  {
    rot_eps_ = r;
  }

  inline void setTransEps(const double t)
  {
    trans_eps_ = t;
  }

  inline void setMaxStepReadjust(const std::size_t s)
  {
    max_step_readjust_ = s;
  }

  inline void setAlpha(const double a)
  {
    alpha_ = a;
  }

  inline void setTransform(const cslibs_math_3d::Transform3d &t)
  {
    transform_ = t;
  }


private:
  double                      resolution_;       /// resolution_
  std::size_t                 max_iterations_;   /// maximum iterations
  double                      rot_eps_;          ///
  double                      trans_eps_;        ///
  std::size_t                 max_step_readjust_;
  double                      alpha_;            /// step correction
  cslibs_math_3d::Transform3d transform_;        /// initial transform_


};
}
}

#endif // CSLIBS_NDT_PARAMS_HPP
