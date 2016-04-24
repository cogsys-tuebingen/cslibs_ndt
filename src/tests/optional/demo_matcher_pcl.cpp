/// PROJECT
#include <ndt/data/pointcloud.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/visualization/visualize.hpp>

#include <pcl/registration/ndt.h>

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

typedef ndt::matching::MultiGridMatcher2D  MatcherType;

int main(int argc, char *argv[])
{
    std::vector<MatcherType::PointType> points_src;
    /// generate horizontal lines
    std::vector<double> xs;
    linspace(-10.0, -1.0, 0.1, xs);
    for(double &e : xs) {
        points_src.push_back(MatcherType::PointType(e, 1.0));
        points_src.push_back(MatcherType::PointType(e, -1.0));
    }
    /// generate vertial lines
    std::vector<double> ys;
    linspace(-10.0, 10.0, 0.1, ys);
    for(double &e : ys) {
        points_src.push_back(MatcherType::PointType(1.5, e));
        if(e < -1.0 || e > 1.0)
            points_src.push_back(MatcherType::PointType(-1.0, e));
    }

    MatcherType::RotationType    rotation       = MatcherType::RotationType(0.2);
    MatcherType::TranslationType trans          = MatcherType::TranslationType(0.2, 0.0);
    MatcherType::TransformType   transformation = trans * rotation;

    pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_points_src(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_points_dst(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ> result;
    std::vector<MatcherType::PointType> points_dst;
    for(MatcherType::PointType &p : points_src) {
        pcl::PointXYZ p_pcl;
        p_pcl.x = p(0);
        p_pcl.y = p(1);
        pcl_points_src->push_back(p_pcl);
        MatcherType::PointType p_transformed = transformation * p;
        p_pcl.x = p_transformed(0);
        p_pcl.y = p_transformed(1);
        pcl_points_dst->push_back(p_pcl);
        points_dst.push_back(p_transformed);
    }

    pcl::NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> ndt;
    ndt.setTransformationEpsilon (0.01);
    ndt.setMaximumIterations (400);
    ndt.setInputSource (pcl_points_dst);
    ndt.setInputTarget (pcl_points_src);
    ndt.align (result);

    /// show the point set
    ndt::MultiGrid2DType::SizeType   size = {20, 20};
    ndt::MultiGrid2DType::ResolutionType resolution = {1.0, 1.0};
    cv::Mat display = cv::Mat(800, 800, CV_8UC3, cv::Scalar());
    ndt::renderPoints(points_src, size, resolution, display, cv::Scalar(255), false, 0.5);
    ndt::renderPoints(points_dst, size, resolution, display, cv::Scalar(0,255), false, 0.5);
    while(true) {
        cv::imshow("display", display);
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }
    cv::flip(display, display, 0);

    points_dst.clear();
    for(pcl::PointXYZ &p_pcl : result) {
        MatcherType::PointType p(p_pcl.x, p_pcl.y);
        points_dst.push_back(p);
    }

    /// now we can try out the matching
    ndt::renderPoints(points_dst, size, resolution, display, cv::Scalar(0,0,255), false, 0.5);
    while(true) {
        cv::imshow("display", display);
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }

    return 0;
}
