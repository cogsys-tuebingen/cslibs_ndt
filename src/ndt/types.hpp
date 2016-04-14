#ifndef TYPES_HPP
#define TYPES_HPP
#include <eigen3/Eigen/Core>

namespace ndt {
typedef Eigen::Vector2d            Point;
typedef Eigen::Matrix<double,1,2>  PointTransposed;
typedef Eigen::Vector3d            Transform;
typedef Eigen::Matrix3d            TransformMatrix;
typedef Eigen::Matrix2d            Covariance;
typedef Eigen::Matrix<double,2,3>  Jacobian;
typedef Eigen::Matrix3d            Hessian;
}
#endif // TYPES_HPP
