#ifndef CSLIBS_NDT_3D_HESSIAN_HPP
#define CSLIBS_NDT_3D_HESSIAN_HPP

#include <Eigen/Eigen>

namespace cslibs_ndt_3d {
namespace matching {
class Hessian {
public:
  using hessian_t = std::array<std::array<Eigen::Matrix3d, 3>, 3>;
  using point_t   = Eigen::Vector3d;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  inline Hessian()
  {
    for(std::size_t i = 0 ; i < 3 ; ++i) {
      for(std::size_t j = 0 ; j < 3 ; ++j) {
        data_[i][j] = Eigen::Matrix3d::Zero();
      }
    }
  }

  enum Partial{tx = 0, ty = 1, tz = 2, alpha = 3, beta = 4, gamma = 5, roll = 3, pitch = 4, yaw = 5};

  inline const Eigen::Vector3d get(const Partial pi,
                                   const Partial pj,
                                   const point_t &p) const
  {
    return (pi < 3 || pj < 3) ? Eigen::Vector3d::Zero() : static_cast<Eigen::Vector3d>(data_[pi-3][pj-3] * p);
  }

  inline const Eigen::Vector3d get(const std::size_t pi,
                                   const std::size_t pj,
                                   const point_t &p) const
  {
    assert(pi < 6);
    assert(pj < 6);
    return (pi < 3 || pj < 3) ? Eigen::Vector3d::Zero() : static_cast<Eigen::Vector3d>(data_[pi-3][pj-3] * p);
  }


  inline const hessian_t & angular() const
  {
    return data_;
  }

  inline hessian_t & angular()
  {
    return data_;
  }


  inline static void get(const std::array<double, 3> &angular,
                         Hessian &h)
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

    std::array<std::array<Eigen::Matrix3d, 3>, 3> &data = h.angular();

    data[0][0](0,1) = -sa*sb*cg + sg*ca;
    data[0][0](0,2) = -sa*sg - sb*ca*cg;
    data[0][0](1,1) = -sa*sb*sg - ca*cg;
    data[0][0](1,2) =  sa*cg - sb*sg*ca;
    data[0][0](2,1) = -sa*cb;
    data[0][0](2,2) = -ca*cb;

    data[0][1](0,1) =  ca*cb*cg;
    data[0][1](0,2) = -sa*cb*cg;
    data[0][1](1,1) =  sg*ca*cb ;
    data[0][1](1,2) = -sa*sg*cb;
    data[0][1](2,1) = -sb*ca;
    data[0][1](2,2) =  sa*sb;

    data[0][2](0,1) =  sa*cg - sb*sg*ca;
    data[0][2](0,2) =  sa*sb*sg + ca*cg;
    data[0][2](1,1) =  sa*sg + sb*ca*cg;
    data[0][2](1,2) = -sa*sb*cg + sg*ca;

    data[1][0](0,1) =  ca*cb*cg;
    data[1][0](0,2) = -sa*cb*cg;
    data[1][0](1,1) =  sg*ca*cb;
    data[1][0](1,2) = -sa*sg*cb;
    data[1][0](2,1) = -sb*ca;
    data[1][0](2,2) =  sa*sb;

    data[1][1](0,0) = -cb*cg;
    data[1][1](0,1) = -sa*sb*cg;
    data[1][1](0,2) = -sb*ca*cg;
    data[1][1](1,0) = -sg*cb;
    data[1][1](1,1) = -sa*sb*sg;
    data[1][1](1,2) = -sb*sg*ca;
    data[1][1](2,0) =  sb;
    data[1][1](2,1) = -sa*cb;
    data[1][1](2,2) = -ca*cb;

    data[1][2](0,0) =  sb*sg;
    data[1][2](0,1) = -sa*sg*cb;
    data[1][2](0,2) = -sg*ca*cb;
    data[1][2](1,0) = -sb*cg;
    data[1][2](1,1) =  sa*cb*cg;
    data[1][2](1,2) =  ca*cb*cg;

    data[2][0](0,1) = sa*cg - sb*sg*ca;
    data[2][0](0,2) = sa*sb*sg + ca*cg;
    data[2][0](1,1) = sa*sg + sb*ca*cg;
    data[2][0](1,2) = -sa*sb*cg + sg*ca;

    data[2][1](0,0) =  sb*sg;
    data[2][1](0,1) = -sa*sg*cb;
    data[2][1](0,2) = -sg*ca*cb;
    data[2][1](1,0) = -sb*cg;
    data[2][1](1,1) =  sa*cb*cg;
    data[2][1](1,2) =  ca*cb*cg;

    data[2][2](0,0) = -cb*cg;
    data[2][2](0,1) = -sa*sb*cg + sg*ca;
    data[2][2](0,2) = -sa*sg - sb*ca*cg;
    data[2][2](1,0) = -sg*cb;
    data[2][2](1,1) = -sa*sb*sg - ca*cg;
    data[2][2](1,2) =  sa*cg - sb*sg*ca;
  }


private:
  hessian_t data_;
};


}
}

#endif // CSLIBS_NDT_3D_HESSIAN_HPP
