#ifndef CSLIBS_NDT_3D_JACOBIAN_HPP
#define CSLIBS_NDT_3D_JACOBIAN_HPP

#include <Eigen/Eigen>

namespace cslibs_ndt_3d {
namespace matching {
class EIGEN_ALIGN16 Jacobian {
public:
    using linear_jacobian_t  = std::array<Eigen::Vector3d, 3>;
    using angular_jacobian_t = std::array<Eigen::Matrix3d, 3>;
    using point_t            = Eigen::Vector3d;
    using matrix_t           = Eigen::Matrix3d;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    inline Jacobian() :
        linear_data_{{Eigen::Vector3d(1.,0.,0.),
                     Eigen::Vector3d(0.,1.,0.),
                     Eigen::Vector3d(0.,0.,1.)}},
        angular_data_{{Eigen::Matrix3d::Zero(),
                      Eigen::Matrix3d::Zero(),
                      Eigen::Matrix3d::Zero()}},
        angular_transposed_data_{{Eigen::Matrix3d::Zero(),
                                 Eigen::Matrix3d::Zero(),
                                 Eigen::Matrix3d::Zero()}}
    {
    }

    enum Partial{tx = 0, ty = 1, tz = 2, alpha = 3, beta = 4, gamma = 5, roll = 3, pitch = 4, yaw = 5};

    inline const point_t get(const Partial  pi,
                             const point_t &p) const
    {
        return  pi < 3 ? linear_data_[pi] : angular_data_[pi - 3] * p;
    }

    inline const point_t get(const std::size_t  pi,
                             const point_t &p) const
    {
        assert(pi < 6);
        return  pi < 3 ? linear_data_[pi] : angular_data_[pi - 3] * p;
    }


    inline const matrix_t get(const Partial  pi,
                              const matrix_t &R,
                              const matrix_t &C) const
    {
        return pi < 3 ? matrix_t::Zeros() : angular_transposed_data_[pi - 3] * C * R +
                                            R * C * angular_data_[pi - 3];
    }

    inline const matrix_t get(const std::size_t pi,
                              const matrix_t &R,
                              const matrix_t &C) const
    {
        assert(pi < 6);
        return pi < 3 ? matrix_t::Zeros() : angular_transposed_data_[pi - 3] * C * R +
                                            R * C * angular_data_[pi - 3];
    }

    inline const angular_jacobian_t & angular() const
    {
        return angular_data_;
    }

    inline angular_jacobian_t & angular()
    {
        return angular_data_;
    }

    //  inline static void get(const std::array<double, 3> &angular, /// linear components not required because the derivation is always the same
    inline static void get(const Eigen::Vector3d &angular, /// linear components not required because the derivation is always the same
                           Jacobian &j)                          /// roll pitch yaw / alpha beta gamma
    {
        const double alpha = angular[0];
        const double beta  = angular[1];
        const double gamma = angular[2];

        const double sa = std::sin(alpha);
        const double sb = std::sin(beta);
        const double sg = std::sin(gamma);
        const double ca = std::cos(alpha);
        const double cb = std::cos(beta);
        const double cg = std::cos(gamma);

        std::array<Eigen::Matrix3d, 3> &data = j.angular_data_;
        std::array<Eigen::Matrix3d, 3> &data_transposed = j.angular_transposed_data_;

        data[0](0,1) =  sa*sg  + sb*ca*cg;
        data[0](0,2) = -sa*sb*cg + sg*ca;
        data[0](1,1) = -sa*cg + sb*sg*ca;
        data[0](1,2) = -sa*sb*sg - ca*cg;
        data[0](2,1) =  ca*cb;
        data[0](2,2) = -sa*cb;

        data[1](0,0) = -sb*cg;
        data[1](0,1) =  sa*cb*cg;
        data[1](0,2) =  ca*cb*cg;
        data[1](1,0) = -sb*sg;
        data[1](1,1) =  sa*sg*cb;
        data[1](1,2) =  sg*ca*cb;
        data[1](2,0) = -cb;
        data[1](2,1) = -sa*sb;
        data[1](2,2) = -sb*ca;

        data[2](0,0) = -sg*cb;
        data[2](0,1) = -sa*sb*sg - ca*cg;
        data[2](0,2) =  sa*cg - sb*sg*ca;
        data[2](1,0) =  cb*cg;
        data[2](1,1) =  sa*sb*cg - sg*ca;
        data[2](1,2) =  sa*sg + sb*ca*cg;

        data_transposed[0] = (data[0].transpose()).eval();
        data_transposed[1] = (data[1].transpose()).eval();
        data_transposed[2] = (data[2].transpose()).eval();
    }

private:
    linear_jacobian_t   linear_data_;
    angular_jacobian_t  angular_data_;
    angular_jacobian_t  angular_transposed_data_;
} ;
}
}
#endif // CSLIBS_NDT_3D_JACOBIAN_HPP
