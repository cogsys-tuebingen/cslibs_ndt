/// PROJECT
#include <ndt/data/pointcloud.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/matching/multi_grid_matcher_2D_ls.hpp>
#include <ndt/visualization/multi_grid.hpp>
#include <ndt/visualization/points.hpp>
#include <ndt/math/hausdorff.hpp>
#include <ndt/matching/kdtree_matcher_2D.hpp>

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

typedef ndt::matching::MultiGridMatcher2DLS MultiGridMatcher2D;
typedef ndt::visualization::MultiGrid2D     MultiGrid2D;
typedef ndt::visualization::Point2D         Point2D;

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
        MultiGridMatcher2D::TranslationType trans          = MultiGridMatcher2D::TranslationType(-0.2, 0.0);
        MultiGridMatcher2D::TransformType   transformation = trans * rotation;
        std::vector<MultiGridMatcher2D::PointType> points_dst;
        for(MultiGridMatcher2D::PointType &p : points_src) {
            points_dst.push_back(transformation * p);
        }

        MultiGrid2D::SizeType   size = {20, 20};
        MultiGridMatcher2D::ResolutionType resolution = {1.0, 1.0};
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

        MultiGridMatcher2D::Parameters params;
        params.eps_rot = 1e-6;
        params.eps_trans = 1e-6;
        params.alpha = 1.45;
        params.lambda = MultiGridMatcher2D::LambdaType::Constant(0.1);
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();
        MultiGridMatcher2D multi_matcher(params);
        multi_matcher.match(pointcloud_dst, pointcloud_src, transformation, prior);
        std::chrono::microseconds elapsed =
                std::chrono::duration_cast<std::chrono::microseconds>
                (std::chrono::system_clock::now() - start);

        std::vector<MultiGridMatcher2D::PointType> points_src_corr = points_src;
        for(MultiGridMatcher2D::PointType &p : points_src_corr) {
            p = transformation * p;
        }


        multi_matcher.printDebugInfo();
        std::cout << "elapsed " << elapsed.count() / 1000.0 << " ms" << std::endl;
        std::cout << "hausdorff " << ndt::math::hausdorff<2>(pointcloud_dst, ndt::data::Pointcloud<2>(points_src_corr)) << std::endl;
        std::cout << "hausdorff_frac " << ndt::math::hausdorffFraction<2>(pointcloud_dst, ndt::data::Pointcloud<2>(points_src_corr), 0.1) << std::endl;
        std::cout << "hausdorff_avg " << ndt::math::hausdorffAvg<2>(pointcloud_dst, ndt::data::Pointcloud<2>(points_src_corr)) << std::endl;
        std::cout << "hausdorff_mpe " << ndt::math::hausdorffMPE<2>(pointcloud_dst, ndt::data::Pointcloud<2>(points_src_corr)) << std::endl;


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
