#include "ndt_map_loader.h"

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "NDTMapLoader");

    cslibs_ndt_2d::NDTMapLoader ml;
    ml.run<double>() || ml.run<float>();

    return 0;
}
