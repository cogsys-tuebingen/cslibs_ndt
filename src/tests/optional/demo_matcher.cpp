/// PROJECT
#include <ndt/data/pointcloud.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/matching/multi_matcher.hpp>
#include <ndt/visualization/multi_grid.hpp>
#include <ndt/visualization/points.hpp>

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

typedef ndt::matching::MultiGridMatcher2D   MultiGridMatcher2D;
typedef ndt::visualization::MultiGrid2D     MultiGrid2D;
typedef ndt::visualization::Point2D         Point2D;
typedef ndt::matching::MultiMatcher<MultiGridMatcher2D> MultiMatcherType;

int main(int argc, char *argv[])
{
    std::vector<MultiGridMatcher2D::PointType> points_src;
    /// generate horizontal lines
    std::vector<double> xs;
    linspace(-10.0, -1.0, 0.1, xs);
    for(double &e : xs) {
        points_src.push_back(MultiGridMatcher2D::PointType(e, 1.0));
        points_src.push_back(MultiGridMatcher2D::PointType(e, -1.0));
    }
    /// generate vertial lines
    std::vector<double> ys;
    linspace(-10.0, 10.0, 0.1, ys);
    for(double &e : ys) {
        points_src.push_back(MultiGridMatcher2D::PointType(1.5, e));
        if(e < -1.0 || e > 1.0)
            points_src.push_back(MultiGridMatcher2D::PointType(-1.0, e));
    }
    MultiGridMatcher2D::TransformType prior = MultiGridMatcher2D::TransformType::Identity();
    for(std::size_t i = 0 ; i < 10 ; ++i) {
        /// generate a second points test which is transformed
        MultiGridMatcher2D::RotationType    rotation       = MultiGridMatcher2D::RotationType(0.1 * i);
        MultiGridMatcher2D::TranslationType trans          = MultiGridMatcher2D::TranslationType(-0.0, 0.0);
        MultiGridMatcher2D::TransformType   transformation = trans * rotation;
        std::vector<MultiGridMatcher2D::PointType> points_dst;
        for(MultiGridMatcher2D::PointType &p : points_src) {
            points_dst.push_back(transformation * p);
        }

        MultiGridMatcher2D::SizeType   size = {10, 10};
        MultiGridMatcher2D::ResolutionType resolution = {5.0, 5.0};
        ndt::data::Pointcloud<2> pointcloud_src(points_src);
        ndt::data::Pointcloud<2> pointcloud_dst(points_dst);

        /// show the point set
        cv::Mat display = cv::Mat(800, 800, CV_8UC3, cv::Scalar());
        ndt::visualization::renderPoints(points_src,
                                         size,
                                         resolution,
                                         display,
                                         cv::Scalar(255),
                                         false, 0.5);
        ndt::visualization::renderPoints(points_dst,
                                         size,
                                         resolution,
                                         display,
                                         cv::Scalar(0,255),
                                         false, 0.5);
        while(true) {
            cv::imshow("display", display);
            int key = cv::waitKey(0) & 0xFF;
            if(key == 27)
                break;
        }

        /// now we can try out the matching

        MultiMatcherType::ParameterSet param_set;
        MultiGridMatcher2D::Parameters params;
        params.max_iterations = 30;
        params.eps_trans = 1e-3;
        params.eps_rot = 1e-3;
        params.lambda(0) = 3;
        params.lambda(1) = 3;
        params.lambda(2) = 1;
        params.resolution = resolution;
        params.max_step_corrections = 10;
        param_set.push_back(params);
        params.resolution[0] *= 0.5;
        params.resolution[1] *= 0.5;
        params.lambda = params.lambda * 0.5;
        param_set.push_back(params);
        params.resolution[0] *= 0.5;
        params.resolution[1] *= 0.5;
        params.lambda = params.lambda * 0.5;
        param_set.push_back(params);
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();
        MultiMatcherType multi_matcher(param_set);
        multi_matcher.match(pointcloud_dst, pointcloud_src, transformation, prior);
        std::chrono::microseconds elapsed =
                std::chrono::duration_cast<std::chrono::microseconds>
                (std::chrono::system_clock::now() - start);
        std::cout << "elapsed " << elapsed.count() / 1000.0 << " ms" << std::endl;

        std::vector<MultiGridMatcher2D::PointType> points_src_corr = points_src;
        for(MultiGridMatcher2D::PointType &p : points_src_corr) {
            p = transformation * p;
        }
        ndt::visualization::renderPoints(points_src_corr,
                                         size,
                                         resolution,
                                         display,
                                         cv::Scalar(0,0,255),
                                         false, 0.5);

        MultiGrid2D grid(size, resolution, pointcloud_dst.min);
        grid.add(points_dst);
        cv::Mat display_distribution(500,500,CV_8UC3, cv::Scalar());
        ndt::visualization::renderMultiGrid(grid, Point2D(-10,-10), Point2D(10,10), display_distribution);

        while(true) {
            cv::imshow("display", display);
            cv::imshow("display_distribution", display_distribution);
            int key = cv::waitKey(0) & 0xFF;
            if(key == 27)
                break;
        }

        prior = transformation;

    }
    return 0;
}
