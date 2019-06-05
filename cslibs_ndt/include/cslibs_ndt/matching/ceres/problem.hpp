#ifndef CSLIBS_NDT_MATCHING_CERES_PROBLEM_HPP
#define CSLIBS_NDT_MATCHING_CERES_PROBLEM_HPP

#include <cslibs_ndt/matching/ceres/translation/translation_cost_functor_2d.hpp>
#include <cslibs_ndt/matching/ceres/translation/translation_cost_functor_3d.hpp>

#include <cslibs_ndt/matching/ceres/rotation/rotation_cost_functor_2d.hpp>
#include <cslibs_ndt/matching/ceres/rotation/rotation_cost_functor_3d_quaternion.hpp>
#include <cslibs_ndt/matching/ceres/rotation/rotation_cost_functor_3d_rpy.hpp>

#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_2d.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_3d_quaternion.hpp>
#include <cslibs_ndt/matching/ceres/map/scan_match_cost_functor_3d_rpy.hpp>

#include <cslibs_ndt/matching/ceres/local_parameterization.hpp>

#include <type_traits>
#include <ceres/problem.h>
#include <ceres/local_parameterization.h>

namespace cslibs_ndt {
namespace matching {
namespace ceres {

template <typename ndt_t, Flag flag_t = Flag::DIRECT, typename ... args_t>
inline void Problem2d(const double& translation_weight, const double& rotation_weight, const double& map_weight,
                      const cslibs_math::linear::Vector<double,2>& translation, const double& rotation,
                      double* ceres_translation, double* ceres_rotation,
                      ::ceres::Problem& problem,
                      const args_t &...args)
{
    problem.AddParameterBlock(ceres_translation, 2, nullptr); //TODO
    problem.AddParameterBlock(ceres_rotation, 1, nullptr);    //TODO

    if (map_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::ScanMatchCostFunctor2dCreator<ndt_t, flag_t>::
                  CreateAutoDiffCostFunction(map_weight, args...),
            nullptr,
            ceres_translation,
            ceres_rotation);
    }
    if (translation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::TranslationCostFunctor2d::
                  CreateAutoDiffCostFunction(translation_weight, translation),
            nullptr,
            ceres_translation);
    }
    if (rotation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::RotationCostFunctor2d::
                  CreateAutoDiffCostFunction(rotation_weight, rotation),
            nullptr,
            ceres_rotation);
    }
}

template <typename ndt_t, Flag flag_t = Flag::DIRECT, typename ... args_t>
inline void Problem3dQuaternion(const double& translation_weight, const double& rotation_weight, const double& map_weight,
                                const cslibs_math::linear::Vector<double,3>& translation, const cslibs_math_3d::Quaterniond& rotation,
                                double* ceres_translation, double* ceres_rotation,
                                ::ceres::Problem& problem,
                                const bool only_yaw,
                                const args_t &...args)
{
    problem.AddParameterBlock(ceres_translation, 3, only_yaw ? new ::ceres::SubsetParameterization(3, { 2 }) :
                                                               nullptr);
    problem.AddParameterBlock(ceres_rotation, 4, only_yaw ? YawOnlyQuaternionPlus::CreateAutoDiff() :
                                                            new ::ceres::QuaternionParameterization());

    if (map_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::ScanMatchCostFunctor3dQuaternionCreator<ndt_t, flag_t>::
                  CreateAutoDiffCostFunction(map_weight, args...),
            nullptr,
            ceres_translation,
            ceres_rotation);
    }
    if (translation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::TranslationCostFunctor3d::
                  CreateAutoDiffCostFunction(translation_weight, translation),
            nullptr,
            ceres_translation);
    }
    if (rotation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::RotationCostFunctor3dQuaternion::
                  CreateAutoDiffCostFunction(rotation_weight, rotation),
            nullptr,
            ceres_rotation);
    }
}

template <typename ndt_t, Flag flag_t = Flag::DIRECT, typename ... args_t>
inline void Problem3dRPY(const double& translation_weight, const double& rotation_weight, const double& map_weight,
                         const cslibs_math::linear::Vector<double,3>& translation, const cslibs_math::linear::Vector<double,3>& rotation,
                         double* ceres_translation, double* ceres_rotation,
                         ::ceres::Problem& problem,
                         const args_t &...args)
{
    problem.AddParameterBlock(ceres_translation, 3, nullptr); //TODO
    problem.AddParameterBlock(ceres_rotation, 3, nullptr);    //TODO

    if (map_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::ScanMatchCostFunctor3dRPYCreator<ndt_t, flag_t>::
                  CreateAutoDiffCostFunction(map_weight, args...),
            nullptr,
            ceres_translation,
            ceres_rotation);
    }
    if (translation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::TranslationCostFunctor3d::
                  CreateAutoDiffCostFunction(translation_weight, translation),
            nullptr,
            ceres_translation);
    }
    if (rotation_weight != 0.0) {
      problem.AddResidualBlock(
            cslibs_ndt::matching::ceres::RotationCostFunctor3dRPY::
                  CreateAutoDiffCostFunction(rotation_weight, rotation),
            nullptr,
            ceres_rotation);
    }
}

}
}
}

#endif // CSLIBS_NDT_MATCHING_CERES_PROBLEM_HPP
