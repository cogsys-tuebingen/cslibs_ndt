#include <ndt/data/laserscan.hpp>
#include <ndt/visualization/points.hpp>
#include <ndt/visualization/kdtree.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/matching/kdtree_matcher_2D.hpp>
#include <ndt/visualization/multi_grid.hpp>

#include <string>

typedef ndt::matching::MultiGridMatcher2D      MatcherType;
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


//    src = dst;
//    MatcherType::RotationType    rotation       = MatcherType::RotationType(0.1);
//    MatcherType::TranslationType trans          = MatcherType::TranslationType(-0.5, 0.0);
//    MatcherType::TransformType   transformation = trans * rotation;
//    for(MatcherType::PointType &p : src.points_data) {
//        p = transformation * p;
//    }

//    srand(0);
//    double rr  = 0.05;
//    for(std::size_t i = 0 ; i < src.size ; ++i) {
//        src.points_data[i](0) += rr - (std::rand() % 1000) / 1000.0 * (rr * 2);
//        src.points_data[i](1) += rr - (std::rand() % 1000) / 1000.0 * (rr * 2);
////        dst.points_data[i] = rot * dst.points_data[i];
////        src.points_data[i] = rot * src.points_data[i];
//    }

    /// now let us render the input, so we can get some overview ;)
    ndt::visualization::Size2D   size = {30, 30};
    ndt::visualization::Resolution2D resolution = {2.0, 2.0};
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

    MatcherType::Parameters params;
    params.max_iterations = 4000;
    params.eps_trans = 1e-10;
    params.eps_rot = 1e-10;
    params.lambda = 2;
    params.resolution = resolution;
    params.max_step_corrections = 10;
    MatcherType matcher(params);
    MatcherType::TransformType  transform = MatcherType::TransformType::Identity();
    double score = matcher.match(dst, src, transform);
    matcher.printDebugInfo();
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
