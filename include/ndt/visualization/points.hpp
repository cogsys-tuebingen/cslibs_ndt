#ifndef POINTS_HPP
#define POINTS_HPP

/// SYSTEM
#include <opencv2/opencv.hpp>

/// PROJECT
#include <ndt/data/pointcloud.hpp>


namespace ndt {
namespace visualization {

typedef data::Pointcloud<2>::PointType  Point2D;
typedef std::array<std::size_t, 2>      Size2D;
typedef std::array<double, 2>           Resolution2D;

void renderPoints(const std::vector<Point2D> &points,
                  const Size2D               &grid_dimension,
                  const Resolution2D         &resolution,
                  cv::Mat                    &dst,
                  const cv::Scalar           &color = cv::Scalar(255),
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


    for(const Point2D &p : points) {
        cv::Point cvp(p(0) * scale_x + dst.cols / 2,
                      dst.rows - 1 - (p(1) * scale_y + dst.rows / 2));
        cv::circle(dst, cvp, 3, color, CV_FILLED, CV_AA);
    }
}
}
}

#endif // POINTS_HPP
