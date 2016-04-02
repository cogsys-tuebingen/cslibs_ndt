#ifndef VISUALIZE_HPP
#define VISUALIZE_HPP
/// PROJECT
#include "../ndt/types.hpp"

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
}

#endif // VISUALIZE_HPP
