#ifndef MATCHER_HPP
#define MATCHER_HPP
#include <memory>
#include "multi_grid.hpp"
#include "../data/pointcloud.hpp"
#include <eigen3/Eigen/Geometry>
#include <stdexcept>

namespace ndt {
template<std::size_t Dim>
class NDTMatcher {
public:
    typedef std::shared_ptr<NDTMatcher>                  Ptr;
    typedef data::Pointcloud<Dim>                        PointCloudType;
    typedef typename PointCloudType::PointType           PointType;
    typedef Eigen::Transform<double, Dim, Eigen::Affine> TransformType;
    typedef NDTMultiGrid<Dim>                            NDTGridType;
    typedef typename NDTMultiGrid<Dim>::Resolution       ResolutionType;
    typedef NDTMatcher<Dim>                              BaseClass;
    typedef typename NDTMultiGrid<Dim>::Matrix           CovarianceMatrixType;
    typedef Eigen::Matrix<double, 2 * Dim - 1, 1>        GradientType;
    typedef typename NDTMultiGrid<Dim>::Distribution     DistributionType;
    typedef typename NDTMultiGrid<Dim>::Distributions    DistributionsType;

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
        typename NDTGridType::Size size;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            if(range(i) <= 0.0)
                throw std::range_error("Point cloud boundaries are not set properly!");
            size[i] = floor(range(i) / resolution[i] + 0.5);
        }

        grid.reset(new NDTGridType(size, resolution, _src.min));
        grid->add(_src);

        return doMatch(_dst, _transformation, _max_iterations, _eps, _eps_rad);
    }

protected:
    typename NDTGridType::Ptr grid;
    ResolutionType            resolution;

    virtual double doMatch(const PointCloudType &_dst,
                           TransformType        &_transformation,
                           const std::size_t     _max_iterations = 100,
                           const double          _eps = 1e-3,
                           const double          _eps_rad = 1e-6)
    {
        /// 2. Initialize the parameter estimate.
        /// 3. For each sample fo the second scan: Map the reconstructed 2D point into the
        ///    the coordinate frame of the first scan according to the parameters.
        /// 4. Determine the corresponding normal distributions for each mapped point.
        /// 5. The score for the parameters is determined by evaluation the distribution
        ///    for each point and summing the result.
        /// 6. Calculate the new parameter estimate by traying to optimize the score.
        ///    This is done by performin one step of Newton's Algorithm.
        /// 7. Go to 3. until a convergence criterion is met.
        return 0.0;
    }

};
}
#endif // MATCHER_HPP
