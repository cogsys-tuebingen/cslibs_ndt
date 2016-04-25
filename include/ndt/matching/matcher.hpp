#ifndef MATCHER_HPP
#define MATCHER_HPP
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>

#include <ndt/data/pointcloud.hpp>
#include <ndt/grid/multi_grid.hpp>

namespace ndt {
namespace matching {
template<std::size_t Dim>
class NDTMatcher {
public:
    typedef std::shared_ptr<NDTMatcher>                  Ptr;
    typedef data::Pointcloud<Dim>                        PointCloudType;
    typedef typename PointCloudType::PointType           PointType;

    typedef grid::MultiGrid<Dim>                GridType;
    typedef typename GridType::ResolutionType       ResolutionType;
    typedef typename GridType::DistributionType     DistributionType;
    typedef typename GridType::DistributionSetType    DistributionsType;
    typedef typename GridType::CovarianceMatrixType           CovarianceMatrixType;

    typedef NDTMatcher<Dim>                              BaseClass;

    typedef Eigen::Matrix<double, 2 * Dim - 1, 1>        GradientType;
    typedef Eigen::Transform<double, Dim, Eigen::Affine> TransformType;

    NDTMatcher(const ResolutionType &_resolution) :
        grid(nullptr),
        resolution(_resolution)
    {
    }

    inline double match(const PointCloudType &_src,
                        const PointCloudType &_dst,
                        TransformType        &_transformation,
                        const std::size_t     _max_iterations = 100,
                        const double          _eps = 1e-3,
                        const double          _eps_rad = 1e-6)
    {
        /// range check & grid size calculation
        /// 1. Build the Normal Distribution Transform of the first scan.
        PointType range = _src.range();
        typename GridType::SizeType size;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            if(range(i) <= 0.0)
                throw std::range_error("Point cloud boundaries are not set properly!");
            size[i] = floor(range(i) / resolution[i] + 0.5);
        }

        grid.reset(new GridType(size, resolution, _src.min));
        grid->add(_src);
        src = _src;

        return doMatch(_dst, _transformation, _max_iterations, _eps, _eps_rad);
    }

protected:
    typename GridType::Ptr    grid;
    PointCloudType            src;
    ResolutionType            resolution;

    virtual double doMatch(const PointCloudType &_dst,
                           TransformType        &_transformation,
                           const std::size_t     _max_iterations = 100,
                           const double          _eps = 1e-3,
                           const double          _eps_rad = 1e-6) = 0;

};
}
}
#endif // MATCHER_HPP
