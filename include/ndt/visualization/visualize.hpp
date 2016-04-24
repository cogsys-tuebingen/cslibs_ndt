#ifndef VISUALIZE_HPP
#define VISUALIZE_HPP
/// PROJECT
#include <ndt/grid/multi_grid.hpp>

/// SYSTEM
#include <opencv2/opencv.hpp>


namespace ndt {

typedef grid::MultiGrid<2> MultiGrid2DType;
typedef grid::Grid<2>      Grid2DType;

void renderPoints(const std::vector<MultiGrid2DType::Point> &points,
                  const MultiGrid2DType::Size               &grid_dimension,
                  const MultiGrid2DType::Resolution         &resolution,
                  cv::Mat &dst,
                  const cv::Scalar &color = cv::Scalar(255),
                  const bool render_grid = true,
                  const double visibility = 0.5)
{
    double scale_x = dst.cols / (grid_dimension[0] * resolution[0]);
    double scale_y = dst.rows / (grid_dimension[1] * resolution[1]);

    if(render_grid) {
        uchar grid_color = visibility * 255;
        for(std::size_t i = 0 ; i < grid_dimension[0] ; ++i) {
            cv::Point a(i * dst.cols / grid_dimension[0], 0);
            cv::Point b(i * dst.cols / grid_dimension[0], dst.rows - 1);
            cv::line(dst, a, b, cv::Scalar(grid_color,grid_color,grid_color), 1, CV_AA);
        }

        for(std::size_t i = 0 ; i < grid_dimension[1] ; ++i) {
            cv::Point a(0, i * dst.rows / grid_dimension[1]);
            cv::Point b(dst.cols - 1, i * dst.rows / grid_dimension[1]);
            cv::line(dst, a, b, cv::Scalar(grid_color,grid_color,grid_color), 1, CV_AA);

        }
    }


    for(const MultiGrid2DType::Point &p : points) {
        cv::Point cvp(p(0) * scale_x + dst.cols / 2,
                      dst.rows - 1 - (p(1) * scale_y + dst.rows / 2));
        cv::circle(dst, cvp, 3, color, CV_FILLED, CV_AA);
    }
}

void renderNDTGrid(MultiGrid2DType              &grid,
                   const MultiGrid2DType::Point &min,
                   const MultiGrid2DType::Point &max,
                   cv::Mat &dst)
{
    double scale_x = fabs((max - min)(0)) / dst.cols;
    double scale_y = fabs((max - min)(1)) / dst.rows;

    cv::Mat samples(dst.rows, dst.cols, CV_64FC1, cv::Scalar());
    double max_value = std::numeric_limits<double>::min();
#pragma omp parallel for reduction(max : max_value)
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            MultiGrid2DType::Point p = min + MultiGrid2DType::Point(scale_x * j, scale_y * i);
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

#endif // VISUALIZE_HPP
