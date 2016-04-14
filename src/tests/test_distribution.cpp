/// SYSTEM
#include <iostream>
#include <chrono>

/// PROJECT
#include "../ndt/rolling_distribution.hpp"
#include "../math/distribution.hpp"

using namespace ndt;

int main(int argc, char *argv[])
{
    Point points[] = { Point(0,0),
                       Point(2,0),
                       Point(2,4),
                       Point(4,4),
                       Point(-8,8)};

    ndt::math::Distribution<2> roll2;
    RollingDistribution roll;
    for(std::size_t i = 0 ; i < 5 ; ++i) {
        roll.add(points[i]);
        roll2.add(points[i]);
    }

    Point mean;
    Covariance cov;
    roll.mean(mean);
    roll.covariance(cov);

    std::cout << mean << std::endl;
    std::cout << cov << std::endl;
    std::cout << "-------" << std::endl;
    std::cout << roll2.getMean() << std::endl;
    std::cout << roll2.getCovariance() << std::endl;


    /// insertion test
    std::chrono::time_point<std::chrono::system_clock> start =
            std::chrono::system_clock::now();
    for(int i = 0 ; i < 5000 ; ++i)
        roll.add(points[0]);
    std::chrono::duration<double> elapsed =
            std::chrono::system_clock::now() - start;
    std::cout << "elapsed " << elapsed.count() << "s" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000.0 << "ms" << std::endl;
    std::cout << "elapsed " << elapsed.count() * 1000000.0 << "µs" << std::endl;

   start = std::chrono::system_clock::now();
    for(int i = 0 ; i < 5000 ; ++i)
        roll2.add(points[0]);
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





    return 0;
}
