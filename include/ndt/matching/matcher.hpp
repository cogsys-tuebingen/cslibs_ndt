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


    struct Parameters {
        Parameters(const ResolutionType &_resolution = {1.0, 1.0},
                   const double          _eps_rot = 1e-6,
                   const double          _eps_trans = 1e-3,
                   const std::size_t     _max_iterations = 100) :
            resolution(_resolution),
            eps_rot(_eps_rot),
            eps_trans(_eps_trans),
            max_iterations(_max_iterations)
        {
        }

        ResolutionType resolution;
        double         eps_rot;
        double         eps_trans;
        std::size_t    max_iterations;
    };

    NDTMatcher(const Parameters &params = Parameters()) :
        grid(nullptr)
    {
    }

    inline double match(const PointCloudType &_src,
                        const PointCloudType &_dst,
                        TransformType        &_transformation)
    {
        /// range check & grid size calculation
        /// 1. Build the Normal Distribution Transform of the first scan.
        PointType range = _src.range();
        typename GridType::SizeType size;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            if(range(i) <= 0.0)
                throw std::range_error("Point cloud boundaries are not set properly!");
            size[i] = floor(range(i) / params.resolution[i] + 0.5);
        }

        grid.reset(new GridType(size, params.resolution, _src.min));
        grid->add(_src);

        return doMatch(_dst, _transformation);
    }

protected:
    typename GridType::Ptr    grid;
    Parameters                params;

    virtual double doMatch(const PointCloudType &_dst,
                           TransformType        &_transformation) = 0;

    inline bool epsTrans(const double a,
                         const double b) const
    {
        return fabs(a - b) < params.eps_trans;
    }

    inline bool epsRot(const double a,
                       const double b) const
    {
        return fabs(a - b) < params.eps_rot;
    }

};
}
}
#endif // MATCHER_HPP
