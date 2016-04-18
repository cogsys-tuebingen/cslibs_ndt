#ifndef MATCHER2D_HPP
#define MATCHER2D_HPP

#include <memory>
#include "multi_grid.hpp"
#include "../data/pointcloud.hpp"
#include <eigen3/Eigen/Geometry>
#include <stdexcept>
#include "matcher.hpp"

namespace ndt {
class NDTMatcher2D : public NDTMatcher<2> {
public:
    typedef Eigen::Matrix<double,2,3>  Jacobian;
    typedef Eigen::Matrix<double,2,3>  Hessian;
    typedef Eigen::Translation2d       Translation;
    typedef Eigen::Rotation2Dd         Rotation;

    NDTMatcher2D(const ResolutionType &resolution) :
        BaseClass(resolution)
    {
    }

private:
    typename NDTGridType::Ptr grid;
    double                    resolution;

    void doMatch(const PointCloudType &_dst,
                 TransformType        &_transformation,
                 const std::size_t     _max_iterations = 100,
                 const double          _eps = 1e-3) override
    {
        _transformation.setIdentity();
        /// 2. Initialize the parameter estimate.
        double theta = 0.0;
        double tx = 0.0;
        double ty = 0.0;
        Translation trans;
        Rotation    rotation(0.0);
        PointType  *points = new PointType[_dst.size];

        /// variables needed for sampling a point
        PointType            mean;
        CovarianceMatrixType inverse_covariance;
        PointType            q;
        double               score;
        PointType            q_inverse_covariance;
        double               sin_theta, cos_theta;
        PointType            jac;

        /// gradient and stuff
        GradientType gradient = GradientType::Zero();

        bool converged = false;
        while(!converged) {
            rotation = Rotation(theta);
            trans    = Translation(tx, ty);
            _transformation = trans * rotation * _transformation;

            for(std::size_t i = 0 ; i < _dst.size ; ++i) {
                if(_dst.mask[i] == PointCloudType::VALID) {
                    points[i] = _transformation * _dst.points[i];
                    score = grid->sampleNonNormalized(points[i], mean, inverse_covariance, q);
                    /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                    /// if s is 0.0 we do no need to count the sample in
                    if(score > 0.0) {
                        q_inverse_covariance = q.transpose() * inverse_covariance;
                        sincos(theta, &sin_theta, &cos_theta);
                        jac(0) = -q(0) * sin_theta - q(1) * cos_theta;
                        jac(1) =  q(0) * cos_theta - q(1) * sin_theta;

                        /// gradient computation
                        gradient(0) -= score * q_inverse_covariance(0);
                        gradient(1) -= score * q_inverse_covariance(1);
                        gradient(2) -= score * double(q_inverse_covariance.dot(jac));
                        /// hessian computation

                    }
                }
            }


        }

        delete[] points;


        /// 3. For each sample fo the second scan: Map the reconstructed 2D point into the
        ///    the coordinate frame of the first scan according to the parameters.
        /// 4. Determine the corresponding normal distributions for each mapped point.
        /// 5. The score for the parameters is determined by evaluation the distribution
        ///    for each point and summing the result.
        /// 6. Calculate the new parameter estimate by traying to optimize the score.
        ///    This is done by performin one step of Newton's Algorithm.
        /// 7. Go to 3. until a convergence criterion is met.

    }
};
}

#endif // MATCHER2D_HPP
