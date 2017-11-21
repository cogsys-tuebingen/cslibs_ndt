#include <eigen3/Eigen/Geometry>
#include <iostream>

int main(int argc, char *argv[])
{
    Eigen::Transform<double, 2, Eigen::Affine> transform;

    Eigen::Vector2d a(0.5, 0.5);
    Eigen::Vector2d b(1.25, 0.25);
    std::cout << a.dot(b) << std::endl;
    std::cout << a.transpose() * b << std::endl;

    return 0;
}
