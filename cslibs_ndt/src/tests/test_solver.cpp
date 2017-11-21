#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

int main()
{
    Matrix3d A = Matrix3d::Identity();
    Vector3d b = Vector3d(1,2,3);
    std::cout << A.jacobiSvd(ComputeFullU | ComputeFullV).solve(b) << endl;
}
