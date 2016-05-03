#include <ndt/data/laserscan.hpp>
#include <ndt/visualization/multi_grid.hpp>

#include <string>

int main(int argc, char *argv[])
{

    if(argc < 3) {
        std::cerr << "ndt_scan_matcher <src> <dst>" << std::endl;
        return 1;
    }

    std::string path_dst = argv[1];
    std::string path_src = argv[2];

    std::cout << "Matching '" << path_src << "' onto '" << path_dst << "'" << std::endl;

    ndt::data::LaserScan src;
    ndt::data::LaserScan dst;
    src.load(path_src);
    dst.load(path_dst);

    std::cout << "ich auch nicht" << std::endl;



    return 0;
}
