#ifndef VISUALIZE_HPP
#define VISUALIZE_HPP
/// PROJECT
#include "../ndt/multi_grid.hpp"

/// SYSTEM
#include <opencv2/opencv.hpp>


namespace ndt {

typedef NDTMultiGrid<2> NDTMultiGrid2D;

void renderPoints(const std::vector<NDTMultiGrid2D::Point> &points,
                  const NDTMultiGrid2D::Size               &grid_dimension,
                  const NDTMultiGrid2D::Resolution         &resolution,
                  cv::Mat &dst,
                  const bool render_grid = true,
                  const double visibility = 0.5,
                  const bool flip = true)
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


    for(const NDTMultiGrid2D::Point &p : points) {
        cv::Point cvp(p(0) * scale_x + dst.cols / 2,
                      p(1) * scale_y + dst.rows / 2);
        cv::circle(dst, cvp, 3, cv::Scalar(255), CV_FILLED, CV_AA);
    }

    if(flip)
        cv::flip(dst, dst, 0);
}

void renderNDTGrid(NDTMultiGrid2D              &grid,
                   const NDTMultiGrid2D::Point &min,
                   const NDTMultiGrid2D::Point &max,
                   cv::Mat &dst)
{
    double scale_x = fabs((max - min)(0)) / dst.cols;
    double scale_y = fabs((max - min)(1)) / dst.rows;

    cv::Mat samples(dst.rows, dst.cols, CV_64FC1, cv::Scalar());
    double max_value = std::numeric_limits<double>::min();
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            NDTMultiGrid2D::Point p = min + NDTMultiGrid2D::Point(scale_x * j, scale_y * i);
            double value = grid.sampleNonNormalized(p);
            if(value > max_value) {
                max_value = value;
            }
            samples.at<double>(i,j) = value;
        }
    }

    cv::imshow("samples", samples);

    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            cv::Vec3b &pix = dst.at<cv::Vec3b>(i,j);
            pix[0] = samples.at<double>(i,j) / max_value * 255;
            pix[1] = samples.at<double>(i,j) / max_value * 255;
            pix[2] = samples.at<double>(i,j) / max_value * 255;
        }
    }
}


//void renderNDTGrid(NDTMultiGrid &multi_grid,
//                   const Point &min,
//                   const Point &max,
//                   const double resolution,
//                   cv::Mat &dst,
//                   const bool flip = true)
//{
//    const double range_x = max(0) - min(0);
//    const double range_y = max(1) - min(1);
//    const cv::Size size(range_x / resolution,
//                        range_y / resolution);
//    const double step_x = range_x / size.width;
//    const double step_y = range_y / size.height;

//    assert(range_x > 0.0);
//    assert(range_y > 0.0);
//    cv::Mat sampled = cv::Mat(size.height, size.width, CV_64FC1, cv::Scalar());
//    Point sample = min;
//    for(int i = 0 ; i < size.height ; ++i) {
//        sample(0)  = min(0);
//        for(int j = 0 ; j < size.width ; ++j) {
//            sample(0) += step_x;
//            sampled.at<double>(i,j) = sqrt(multi_grid.sample(sample));
//        }
//        sample(1) += step_y;
//    }

//    cv::normalize(sampled, sampled, cv::NORM_MINMAX, 0, 1);
//    sampled.convertTo(dst, CV_8UC1, 255.0);
//    cv::cvtColor(dst, dst, CV_GRAY2BGR);
//    const int every_px = multi_grid.resolution() / resolution;
//    for(int i = 0 ; i < size.width ; i += every_px) {
//        for(int j = 0 ; j < size.height; ++j) {
//            cv::Vec3b &c = dst.at<cv::Vec3b>(j,i);
//            c[0] = 255;
//            c[1] = 127;
//            c[2] = 127;
//        }
//    }
//    for(int i = 0 ; i < size.height ; i+=every_px) {
//        for(int j = 0 ; j < size.width ; ++j) {
//            cv::Vec3b &c = dst.at<cv::Vec3b>(i,j);
//            c[0] = 255;
//            c[1] = 127;
//            c[2] = 127;
//        }
//    }
//    if(flip)
//        cv::flip(dst, dst, 0);
//}
}

#endif // VISUALIZE_HPP
