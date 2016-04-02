/// PROJECT
#include "../ndt/ndt_multi_grid.h"

#include "visualize.hpp"

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


int main(int argc, char *argv[])
{
    std::vector<ndt::Point> points;
    /// generate horizontal lines
    std::vector<double> xs;
    linspace(-10.0, -1.0, 0.1, xs);
    for(double &e : xs) {
        points.push_back(ndt::Point(e, 1.0));
        points.push_back(ndt::Point(e, -1.0));
    }
    /// generate vertial lines
    std::vector<double> ys;
    linspace(-10.0, 10.0, 0.1, ys);
    for(double &e : ys) {
        points.push_back(ndt::Point(1.5, e));
        if(e < -1.0 || e > 1.0)
            points.push_back(ndt::Point(-1.0, e));
    }

    cv::Mat display;
    ndt::renderPoints(points, cv::Size(800,800), display);
    cv::imshow("display", display);
    cv::waitKey(0);

    /// no we put the multi grid into play
    ndt::NDTMultiGrid multigrid(20.0, 20.0, 1.0);
    for(ndt::Point &p : points) {
       if(!multigrid.add(p))
           std::cerr << "could not add point" << std::endl;
    }

    cv::Mat display_ndt;
    ndt::renderNDTGrid(multigrid, ndt::Point(-10, -10), ndt::Point(10, 10), 0.01, display_ndt);
    cv::imshow("display_ndt", display_ndt);
    while(true) {
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }
    return 0;
}
