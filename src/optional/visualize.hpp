#ifndef VISUALIZE_HPP
#define VISUALIZE_HPP
/// PROJECT
#include "../ndt/types.hpp"
#include "../ndt/ndt_multi_grid.h"

/// SYSTEM
#include <opencv2/opencv.hpp>


namespace ndt {
void renderPoints(const std::vector<ndt::Point> &points,
                  const cv::Size &size,
                  cv::Mat &dst)
{
    /// find min and max
    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    for(const ndt::Point &p : points) {
        if(p(0) < xmin)
            xmin = p(0);
        if(p(0) > xmax)
            xmax = p(0);
        if(p(1) < ymin)
            ymin = p(1);
        if(p(1) > ymax)
            ymax = p(1);
    }

    double range_x = xmax - xmin;
    double range_y = ymax - ymin;
    double resx = range_x / size.width;
    double resy = range_y / size.height;

    dst = cv::Mat(size.height, size.width, CV_8UC3, cv::Scalar());
    for(const ndt::Point &p : points) {
        cv::Point cvp(p(0) / resx + size.width / 2,
                      p(1) / resy + size.height / 2);
        cv::circle(dst, cvp, 3, cv::Scalar(255), CV_FILLED, CV_AA);
    }

    cv::flip(dst, dst, 0);
}

void renderNDTGrid(NDTMultiGrid &multi_grid,
                   const Point &min,
                   const Point &max,
                   const double resolution,
                   cv::Mat &dst)
{
    const double range_x = max(0) - min(0);
    const double range_y = max(1) - min(1);
    const cv::Size size(range_x / resolution,
                        range_y / resolution);
    const double step_x = range_x / size.width;
    const double step_y = range_y / size.height;

    assert(range_x > 0.0);
    assert(range_y > 0.0);
    cv::Mat sampled = cv::Mat(size.height, size.width, CV_64FC1, cv::Scalar());
    Point sample = min;
    for(int i = 0 ; i < size.height ; ++i) {
        sample(0)  = min(0);
        for(int j = 0 ; j < size.width ; ++j) {
            sample(0) += step_x;
            sampled.at<double>(i,j) = sqrt(multi_grid.sample(sample));
        }
        sample(1) += step_y;
    }

    cv::normalize(sampled, sampled, cv::NORM_L2, 0, 1);
    sampled.convertTo(dst, CV_8UC1, 255.0);
    cv::cvtColor(dst, dst, CV_GRAY2BGR);
    const int every_px = multi_grid.resolution() / resolution;
    for(int i = 0 ; i < size.width ; i += every_px) {
        for(int j = 0 ; j < size.height; ++j) {
            cv::Vec3b &c = dst.at<cv::Vec3b>(j,i);
            c[0] = 255;
            c[1] = 127;
            c[2] = 127;
        }
    }
    for(int i = 0 ; i < size.height ; i+=every_px) {
        for(int j = 0 ; j < size.width ; ++j) {
            cv::Vec3b &c = dst.at<cv::Vec3b>(i,j);
            c[0] = 255;
            c[1] = 127;
            c[2] = 127;
        }
    }



}
}

#endif // VISUALIZE_HPP
