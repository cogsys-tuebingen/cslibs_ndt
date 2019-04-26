#ifndef CSLIBS_NDT_2D_HESSIAN_HPP
#define CSLIBS_NDT_2D_HESSIAN_HPP

#include <Eigen/Eigen>

namespace cslibs_ndt_2d {
namespace matching {
class EIGEN_ALIGN16 Hessian {
public:
    using matrix_t  = Eigen::Matrix2d;
    using hessian_t = std::array<std::array<matrix_t, 1>, 1>;
    using point_t   = Eigen::Vector2d;
    using angular_t = Eigen::Matrix<double, 1, 1>;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline Hessian()
    {
        for(std::size_t i = 0 ; i < 1 ; ++i) {
            for(std::size_t j = 0 ; j < 1 ; ++j) {
                data_[i][j] = matrix_t::Zero();
                data_transposed_[i][j] = matrix_t::Zero();
            }
        }
    }

    enum Partial{tx = 0, ty = 1, gamma = 2, yaw = 2};

    inline const point_t get(const Partial pi,
                             const Partial pj,
                             const point_t &p) const
    {
        return (pi < 2 || pj < 2) ? Eigen::Vector2d::Zero() : static_cast<Eigen::Vector2d>(data_[pi-2][pj-2] * p);
    }

    inline const point_t get(const std::size_t pi,
                             const std::size_t pj,
                             const point_t &p) const
    {
        assert(pi < 3);
        assert(pj < 3);
        return (pi < 2 || pj < 2) ? Eigen::Vector2d::Zero() : static_cast<Eigen::Vector2d>(data_[pi-2][pj-2] * p);
    }

    inline const matrix_t get(const Partial pi,
                              const Partial pj,
                              const matrix_t &C) const
    {
        auto eval = [this, pi, pj, &C]()
        {
            const matrix_t &h = data_[pi-2][pj-2];
            const matrix_t &h_t = data_transposed_[pi-2][pj-2];
            return static_cast<matrix_t>((h_t * C * rotation_).eval() + (rotation_transposed_ * C * h).eval());
        };

        return (pi < 2 || pj < 2) ? matrix_t::Zero() : eval();
    }

    inline const matrix_t get(const std::size_t pi,
                              const std::size_t pj,
                              const matrix_t &C) const
    {
        assert(pi < 3);
        assert(pj < 3);
        auto eval = [this, pi, pj, &C]()
        {
            const matrix_t &h   = data_[pi-2][pj-2];
            const matrix_t &h_t = data_transposed_[pi-2][pj-2];
            return static_cast<matrix_t>((h_t * C * rotation_).eval() + (rotation_transposed_* C * h).eval());
        };

        return (pi < 2 || pj < 2) ? matrix_t::Zero() : eval();
    }


    inline const hessian_t & angular() const
    {
        return data_;
    }

    inline hessian_t & angular()
    {
        return data_;
    }

    inline const matrix_t & rotation() const
    {
        return rotation_;
    }

    inline matrix_t & rotation()
    {
        return rotation_;
    }

    inline static void get(const angular_t &angular,
                           Hessian &h)
    {
        const double &gamma = angular.value();

        const double s = std::sin(gamma);
        const double c = std::cos(gamma);

        hessian_t &data             = h.data_;
        hessian_t &data_transposed  = h.data_transposed_;
        matrix_t  &R                = h.rotation_;
        matrix_t  &R_transposed     = h.rotation_transposed_;

        data[0][0](0,0) = -c;  // dJ(0)/(dyaw,dtx)
        data[0][0](0,1) = s;   // dJ(0)/(dyaw,dty)
        data[0][0](1,0) = -s;  // dJ(1)/(dyaw,dtx)
        data[0][0](1,1) = -c;  // dJ(2)/(dyaw,dty)
        data_transposed[0][0] = data[0][0].transpose();

        R(0,0) = c;
        R(0,1) = -s;
        R(1,0) = s;
        R(1,1) = c;
        R_transposed = R.transpose();
    }


private:
    hessian_t data_;
    hessian_t data_transposed_;
    matrix_t  rotation_;
    matrix_t  rotation_transposed_;
};


}
}

#endif // CSLIBS_NDT_2D_HESSIAN_HPP
