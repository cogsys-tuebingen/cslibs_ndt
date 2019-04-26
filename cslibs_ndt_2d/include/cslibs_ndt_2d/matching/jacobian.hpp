#ifndef CSLIBS_NDT_2D_JACOBIAN_HPP
#define CSLIBS_NDT_2D_JACOBIAN_HPP

#include <eigen3/Eigen/Eigen>

namespace cslibs_ndt_2d {
namespace matching {
class EIGEN_ALIGN16 Jacobian {
public:
    using linear_jacobian_t  = std::array<Eigen::Vector2d, 2>;
    using angular_jacobian_t = std::array<Eigen::Matrix2d, 1>;
    using point_t            = Eigen::Vector2d;
    using matrix_t           = Eigen::Matrix2d;
    using angular_t          = Eigen::Matrix<double, 1, 1>;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline Jacobian() :
        linear_data_{{point_t(1.,0.),
                      point_t(0.,1.)}},
        angular_data_{{matrix_t::Zero()}},
        angular_transposed_data_{{matrix_t::Zero()}},
        linear_rotation_derivative_(matrix_t::Zero()),
        rotation_(matrix_t::Zero())
    {
    }

    enum Partial{tx = 0, ty = 1, gamma = 2, yaw = 2};

    inline const point_t get(const Partial  pi,
                             const point_t &p) const
    {
        return  pi < 2 ? linear_data_[pi] : angular_data_[pi - 2] * p;
    }

    inline const point_t get(const std::size_t  pi,
                             const point_t &p) const
    {
        assert(pi < 3);
        return  pi < 2 ? linear_data_[pi] : angular_data_[pi - 2] * p;
    }


    inline const matrix_t get(const Partial  pi,
                              const matrix_t &C) const
    {
        return pi < 2 ? linear_rotation_derivative_ :
                        (angular_transposed_data_[pi - 2] * C * rotation_).eval() +
                        (rotation_transposed_ * C * angular_data_[pi - 2]).eval();
    }

    inline const matrix_t get(const std::size_t pi,
                              const matrix_t &C) const
    {
        assert(pi < 3);
        return pi < 2 ? linear_rotation_derivative_ :
                        (angular_transposed_data_[pi - 2] * C * rotation_).eval() +
                        (rotation_transposed_ * C * angular_data_[pi - 2]).eval();
    }

    inline const angular_jacobian_t & angular() const
    {
        return angular_data_;
    }

    inline angular_jacobian_t & angular()
    {
        return angular_data_;
    }

    inline const matrix_t& rotation() const
    {
        return rotation_;
    }

    inline static void get(const angular_t &angular,
                           Jacobian &j)
    {
        const double &gamma = angular.value();

        const double s = std::sin(gamma);
        const double c = std::cos(gamma);

        angular_jacobian_t &data            = j.angular_data_;
        angular_jacobian_t &data_transposed = j.angular_transposed_data_;
        matrix_t           &R               = j.rotation_;
        matrix_t           &R_transposed    = j.rotation_transposed_;

        data[0](0,0) = -s;  // dJ(0)/dtx
        data[0](0,1) = -c;  // dJ(0)/dty
        data[0](1,0) = c;   // dJ(1)/dtx
        data[0](1,1) = -s;  // dJ(1)/dty
        data_transposed[0] = (data[0].transpose()).eval();

        R(0,0) = c;
        R(0,1) = -s;
        R(1,0) = s;
        R(1,1) = c;
        R_transposed = R.transpose();
    }

private:
    linear_jacobian_t   linear_data_;
    angular_jacobian_t  angular_data_;
    angular_jacobian_t  angular_transposed_data_;
    matrix_t            linear_rotation_derivative_;
    matrix_t            rotation_;
    matrix_t            rotation_transposed_;
} ;
}
}
#endif // CSLIBS_NDT_2D_JACOBIAN_HPP
