#include "../ndt/matcher2D.hpp"
#include <eigen3/Eigen/Geometry>


int main(int argc, char *argv[])
{
    ndt::NDTMatcher<2> match({1.0, 1.0});

    Eigen::Transform<double, 5, Eigen::Affine> transform;



    return 0;
}
