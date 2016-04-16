#ifndef MATCHER_HPP
#define MATCHER_HPP

#include <memory>
#include "multi_grid.hpp"
#include "../data/pointcloud.hpp"
#include <eigen3/Eigen/Geometry>

namespace ndt {
template<std::size_t Dim>
class NDTMatcher {
public:
    typedef std::shared_ptr<NDTMatcher<Dim>>              Ptr;
    typedef data::Pointcloud<Dim>                         PointCloudType;
    typedef typename PointCloudType::PointType            PointType;
    typedef Eigen::Transform<double, Dim, Eigen::Affine>  TransformType;

    NDTMatcher() :
        grid(nullptr)
    {
    }


private:
    typename NDTMultiGrid<Dim>::Ptr grid;


};


}

#endif // MATCHER_HPP
