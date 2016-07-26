/// PROJECT
#include <ndt/visualization/points.hpp>
#include <ndt/visualization/multi_grid.hpp>
#include <ndt/data/pointcloud.hpp>
#include <ndt/visualization/kdtree.hpp>

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

typedef ndt::visualization::MultiGrid2D  MultiGrid2D;
typedef ndt::visualization::Point2D      Point2D;
typedef ndt::visualization::KDIndex2D    KDIndex2D;
typedef ndt::visualization::KDNodeData2D KDNodeData2D;
typedef ndt::visualization::KDTree2D     KDTree2D;
typedef ndt::visualization::KDNode2D     KDNode2D;

int main(int argc, char *argv[])
{
    std::vector<Point2D> points;
    /// generate horizontal lines
    std::vector<double> xs;
    linspace(-10.0, -1.0, 0.1, xs);
    for(double &e : xs) {
        points.push_back(Point2D(e, 0.0));
        points.push_back(Point2D(e, -2.0));
    }
    /// generate vertial lines
    std::vector<double> ys;
    linspace(-10.0, 10.0, 0.1, ys);
    for(double &e : ys) {
        points.push_back(Point2D(1.5, e));
        if(e < -2.0 || e > 0.0)
            points.push_back(Point2D(-1.0, e));
    }

    ndt::data::Pointcloud<2> pointcloud(points);
    std::cout << "min " << pointcloud.min << std::endl;
    std::cout << "max " << pointcloud.max << std::endl;

    MultiGrid2D::SizeType   size = {20, 20};
    MultiGrid2D::ResolutionType resolution = {1.0, 1.0};

    cv::Mat display = cv::Mat(800, 800, CV_8UC3, cv::Scalar());
    ndt::visualization::renderPoints(points, size, resolution, display);
    cv::imshow("display", display);
    cv::waitKey(0);

    /// no we put the multi grid into play
    MultiGrid2D multigrid(size, resolution, Point2D(-10, -10));
    for(Point2D &p : points) {
        if(!multigrid.add(p))
            std::cerr << "could not add point" << std::endl;
    }
    if(!multigrid.add(pointcloud))
        std::cerr << "some points could not be added" << std::endl;

    ndt::visualization::renderMultiGrid(multigrid, Point2D(-10, -10), Point2D(10,10), display);
    cv::imshow("display", display);
    while(true) {
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }

    KDIndex2D index({2.0, 2.0});
    KDTree2D::Ptr  tree(new KDTree2D);
    for(Point2D &p : points) {
        tree->insert_bulk(index.create(p), p);
    }
    tree->load_bulk();

    ndt::visualization::renderTree(tree, index, Point2D(-10, -10), Point2D(10,10), display);
    cv::imshow("display", display);
    while(true) {
        int key = cv::waitKey(0) & 0xFF;
        if(key == 27)
            break;
    }

    return 0;
}
