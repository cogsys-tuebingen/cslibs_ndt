#include <ndt/data/laserscan.hpp>
#include <ndt/visualization/points.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/visualization/multi_grid.hpp>

#include <string>

typedef ndt::matching::MultiGridMatcher2D   MatcherType;
typedef ndt::visualization::MultiGrid2D     MultiGrid2D;
typedef ndt::visualization::Point2D         Point2D;

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

    MatcherType::RotationType rot(M_PI);

    srand(0);
    double rr  = 0.05;
    for(std::size_t i = 0 ; i < src.size ; ++i) {
//        src.points_data[i](0) += rr - (std::rand() % 1000) / 1000.0 * (rr * 2);
//        src.points_data[i](1) += rr - (std::rand() % 1000) / 1000.0 * (rr * 2);
//        dst.points_data[i] = rot * dst.points_data[i];
//        src.points_data[i] = rot * src.points_data[i];
    }

    /// now let us render the input, so we can get some overview ;)
    ndt::visualization::Size2D   size = {30, 30};
    ndt::visualization::Resolution2D resolution = {1.0, 1.0};
    cv::Mat display(800,800,CV_8UC3,cv::Scalar());
    ndt::visualization::renderPoints(src.points_data,
                                     size,
                                     resolution,
                                     display,
                                     cv::Scalar(255),
                                     false,
                                     0.5);
    ndt::visualization::renderPoints(dst.points_data,
                                     size,
                                     resolution,
                                     display,
                                     cv::Scalar(0, 255),
                                     false,
                                     0.5);

    cv::imshow("display", display);
    while(true) {
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }

    MatcherType matcher;
    MatcherType::TransformType  transform = MatcherType::TransformType::Identity();
    double score = matcher.match(dst, src, transform);
    std::cout << "-------------------------------" << std::endl;
    std::cout << score << std::endl;
    std::cout << transform.translation() << std::endl;
    std::cout << transform.rotation() << std::endl;
    std::cout << "-------------------------------" << std::endl;

    for(auto &p : src.points_data) {
        p = transform * p;
    }

    display = cv::Mat(800,800,CV_8UC3,cv::Scalar());
    ndt::visualization::renderPoints(src.points_data,
                                     size,
                                     resolution,
                                     display,
                                     cv::Scalar(255),
                                     false,
                                     0.5);
    ndt::visualization::renderPoints(dst.points_data,
                                     size,
                                     resolution,
                                     display,
                                     cv::Scalar(0, 255),
                                     false,
                                     0.5);

    cv::imshow("display", display);
    while(true) {
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }


    return 0;
}
