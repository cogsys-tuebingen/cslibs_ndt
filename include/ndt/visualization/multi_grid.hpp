#ifndef VIS_MULTI_GRID_HPP
#define VIS_MULTI_GRID_HPP
/// PROJECT
#include <ndt/grid/multi_grid.hpp>

/// SYSTEM
#include <opencv2/opencv.hpp>
namespace ndt {
namespace visualization {

typedef grid::MultiGrid<2>     MultiGrid2D;
typedef MultiGrid2D::PointType Point2D;

void renderMultiGrid(MultiGrid2D   &grid,
                     const Point2D &min,
                     const Point2D &max,
                     cv::Mat       &dst)
{
    double scale_x = fabs((max - min)(0)) / dst.cols;
    double scale_y = fabs((max - min)(1)) / dst.rows;

    cv::Mat samples(dst.rows, dst.cols, CV_64FC1, cv::Scalar());
    double max_value = std::numeric_limits<double>::lowest();
#pragma omp parallel for
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            MultiGrid2D::PointType p = min + MultiGrid2D::PointType(scale_x * j, scale_y * i);
            double value = grid.sampleNonNormalized(p);
            if(value > max_value) {
                max_value = value;
            }
            samples.at<double>(dst.rows - 1 - i,j) = value;
        }
    }

#pragma omp parallel for
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            cv::Vec3b &pix = dst.at<cv::Vec3b>(i,j);
            pix[0] = samples.at<double>(i,j) / max_value * 255;
            pix[1] = samples.at<double>(i,j) / max_value * 255;
            pix[2] = samples.at<double>(i,j) / max_value * 255;
        }
    }
}
}
}
#endif // VIS_MULTI_GRID_HPP
