/// SYSTEM
#include <iostream>
#include <chrono>

/// PROJECT
#include <ndt/math/distribution.hpp>

using namespace ndt;

void linspace(const double min,
              const double max,
              const double res,
              std::vector<double> &values)
{
    const double        range = max - min;
    const std::size_t   intervals = range / res;
    for(std::size_t i = 0 ; i < intervals ; ++i) {
        values.push_back(min + res * i);
    }
}

typedef Eigen::Vector2d            Point;
typedef Eigen::Matrix<double,1,2>  PointTransposed;
typedef Eigen::Vector3d            Transform;
typedef Eigen::Matrix3d            TransformMatrix;
typedef Eigen::Matrix2d            Covariance;
typedef Eigen::Matrix<double,2,3>  Jacobian;
typedef Eigen::Matrix3d            Hessian;

int main(int argc, char *argv[])
{
    Point points[] = { Point(0,0),
                       Point(2,0),
                       Point(2,4),
                       Point(4,4),
                       Point(-8,8)};

    ndt::math::Distribution<2, false> roll;
    for(std::size_t i = 0 ; i < 5 ; ++i) {
        roll.add(points[i]);
    }
    std::cout << "-------" << std::endl;
    std::cout << roll.getMean() << std::endl;
    std::cout << roll.getCovariance() << std::endl;
    std::cout << "inv_cov" << std::endl;
    std::cout << roll.getInverseCovariance() << std::endl;

    /// insertion test
    std::chrono::time_point<std::chrono::system_clock> start =
            std::chrono::system_clock::now();
    for(int i = 0 ; i < 5000 ; ++i)
        roll.add(points[0]);
    std::chrono::duration<double> elapsed =
            std::chrono::system_clock::now() - start;
    std::cout << "elapsed " << elapsed.count() << "s" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000.0 << "ms" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000000.0 << "ns" << std::endl;

   start = std::chrono::system_clock::now();
    for(int i = 0 ; i < 5000 ; ++i)
        roll.add(points[0]);
    elapsed = std::chrono::system_clock::now() - start;
    std::cout << "elapsed " << elapsed.count() << "s" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000.0 << "ms" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000000.0 << "µs" << std::endl;

    Eigen::Matrix2d test_matrix;
    test_matrix(0,0) = 1;
    test_matrix(0,1) = 2;
    test_matrix(1,0) = 3;
    test_matrix(1,1) = 4;
    for(std::size_t i = 0 ; i < 4 ; ++i)
        std::cout << test_matrix.data()[i] << std::endl;


    std::vector<Point> point_list;
    /// generate horizontal lines
    std::vector<double> xs;
    linspace(-10.0, -1.0, 0.1, xs);
    for(double &e : xs) {
        point_list.push_back(Point(e, 1.0));
        point_list.push_back(Point(e, -1.0));
    }
    /// generate vertial lines
    std::vector<double> ys;
    linspace(-10.0, 10.0, 0.1, ys);
    for(double &e : ys) {
        point_list.push_back(Point(1.5, e));
        if(e < -1.0 || e > 1.0)
            point_list.push_back(Point(-1.0, e));
    }

    roll.reset();
    ndt::math::Distribution<2, false> roll2;
    for(std::size_t i = 0 ; i < 5 ; ++i) {
        roll.add(points[i]);
    }
    for(std::size_t i = 0 ; i < 4 ; ++i) {
        roll2.add(points[i]);
    }
    std::cout << "-------" << std::endl;
    std::cout << roll.getN() << std::endl;
    std::cout << roll.getMean() << std::endl;
    std::cout << roll.getCovariance() << std::endl;
    std::cout << roll.getCovariance().determinant() << std::endl;
    std::cout << roll.sampleNonNoramlized(Point(0.0, 0.0)) << std::endl;

    roll += roll2;
    std::cout << "-------" << std::endl;
    std::cout << roll.getN() << std::endl;
    std::cout << roll.getMean() << std::endl;
    std::cout << roll.getCovariance() << std::endl;
    std::cout << roll.getCovariance().determinant() << std::endl;
    std::cout << roll.sampleNonNoramlized(Point(0.0, 0.0)) << std::endl;


    return 0;
}
