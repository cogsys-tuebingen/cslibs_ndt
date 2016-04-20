
#ifndef MATCHER2D_HPP
#define MATCHER2D_HPP

#include "multi_grid.hpp"
#include "matcher.hpp"
#include "../data/pointcloud.hpp"

#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace ndt {
class NDTMatcher2D : public NDTMatcher<2> {
public:
    typedef Eigen::Matrix<double,3,3>  HessianType;
    typedef Eigen::Translation2d       TranslationType;
    typedef Eigen::Rotation2Dd         RotationType;

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
        /// todo:: initialize parameter estimate
        double phi = 0.0;
        double tx = 0.0;
        double ty = 0.0;
        TranslationType trans;
        RotationType    rotation(0.0);
        PointType  *points = new PointType[_dst.size];

        /// variables needed for sampling a point
        PointType            mean;
        CovarianceMatrixType inverse_covariance;
        PointType            q;
        double               s;
        double               g_dot;
        PointType            q_inverse_covariance;
        double               sin_theta;
        double               cos_theta;
        PointType            jac;
        PointType            hes;


        /// gradient and stuff
        GradientType gradient;
        GradientType delta_p;
        HessianType  hessian;
        double       score;
        double       tx_old;
        double       ty_old;
        double       theta_old;
        std::size_t  iteration = 0;

        bool converged = false;
        while(!converged) {
            rotation = RotationType(phi);
            trans    = TranslationType(tx, ty);
            _transformation = trans * rotation * _transformation;

            gradient = GradientType::Zero();
            hessian  = HessianType::Zero();
            score = 0.0;
            tx_old = tx;
            ty_old = ty;
            theta_old = phi;

            for(std::size_t i = 0 ; i < _dst.size ; ++i) {
                if(_dst.mask[i] == PointCloudType::VALID) {
                    points[i] = _transformation * _dst.points[i];
                    s = grid->sampleNonNormalized(points[i], mean, inverse_covariance, q);
                    /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                    /// if s is 0.0 we do no need to count the sample in
                    if(s > 0.0) {
                        score += s;
                        q_inverse_covariance = q.transpose() * inverse_covariance;
                        sincos(phi, &sin_theta, &cos_theta);
                        jac(0) = -q(0) * sin_theta - q(1) * cos_theta;
                        jac(1) =  q(0) * cos_theta - q(1) * sin_theta;
                        hes(0) = -q(0) * cos_theta + q(1) * sin_theta;
                        hes(1) = -q(0) * sin_theta - q(1) * cos_theta;

                        /// gradient computation
                        g_dot = q_inverse_covariance.dot(jac);
                        gradient(0) -= s * q_inverse_covariance(0);
                        gradient(1) -= s * q_inverse_covariance(1);
                        gradient(2) -= s * g_dot;
                        /// hessian computation
                        hessian(0,0)+=  s
                                     * -q_inverse_covariance(0) * q_inverse_covariance(0)   /// (1)
                                     +  inverse_covariance(0,0);                            /// (3)
                        hessian(1,0)+=  s
                                     * -q_inverse_covariance(1) * q_inverse_covariance(0)   /// (1)
                                     +  inverse_covariance(1,0);                            /// (3)
                        hessian(2,0)+=  s
                                     * -g_dot * inverse_covariance(0)                       /// (1)
                                     +  inverse_covariance.row(0) * jac;                    /// (3)
                        hessian(0,1)+=  s
                                     * -q_inverse_covariance(0) * q_inverse_covariance(1)   /// (1)
                                     +  inverse_covariance(0,1);                            /// (3)
                        hessian(1,1)+=  s
                                     * -q_inverse_covariance(1) * q_inverse_covariance(1)   /// (1)
                                     +  inverse_covariance(1,1);                            /// (3)
                        hessian(2,1)+=  s
                                     * -g_dot * inverse_covariance(1)                       /// (1)
                                     +  inverse_covariance.row(1) * jac;                    /// (3)
                        hessian(0,2)+=  s
                                     * -q_inverse_covariance(0) * g_dot                     /// (1)
                                     +  jac.transpose() * inverse_covariance.col(0);        /// (3)
                        hessian(1,2)+=  s
                                     * -q_inverse_covariance(1) * g_dot                     /// (1)
                                     +  jac.transpose() * inverse_covariance.col(1);        /// (3)
                        hessian(2,2)+=  s
                                     * -g_dot * g_dot                                       /// (1)
                                     +  q_inverse_covariance.dot(hes)                       /// (2)
                                     +  jac.transpose() * inverse_covariance * jac;         /// (3)

                        /// (1) directly computed from Jacobian combined with q^t * InvCov
                        /// (2) only a result for H(2,2)
                        /// (3) [1,0].[[a,b],[c,d]].[[1],[0]] = a with i = j = 0, Jac.col(0)
                        ///     [0,1].[[a,b],[c,d]].[[1],[0]] = c with i = 1, j = 0, Jac.col(1), Jac.col(0)
                        ///     [1,0].[[a,b],[c,d]].[[0],[1]] = b with i = 0, j = 1, Jac.col(0), Jac.col(1)
                        ///     [0,1].[[a,b],[c,d]].[[0],[1]] = d with i = 1, j = 1, Jac.col(1)
                        ///     [1,0].[[a,b],[c,d]].J_T.col(2) = [a,b].J_T.col(2)
                        ///     [0,1].[[a,b],[c,d]].J_T.col(2) = [c,d].J_T.col(2)
                        ///     J_T.col(2).[[a,b],[c,d]].[[1],[0]] = J_T.col(2).[a, c]
                        ///     J_T.col(2).[[a,b],[c,d]].[[0],[1]] = J_T.col(2).[b, d]
                        ///     J_T.col(2).[[a,b],[c,d]].J_T.col(2)
                    }
                }
            }
            /// insert positive definite gurantee here

            /// solve equeation here
            delta_p = GradientType::Zero();
            delta_p = hessian.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(gradient);
            tx  += delta_p(0);
            ty  += delta_p(1);
            phi += delta_p(2);

            /// check for convergence
            if((eps(tx, tx_old, _eps) &&
                    eps(ty, ty_old, _eps) &&
                        eps(phi, theta_old, _eps)) ||
                            iteration > _max_iterations)
                break;
            ++iteration;
        }

        std::cout << "score was " << score << std::endl;

        delete[] points;

    }

    inline bool eps(const double a,
                    const double b,
                    const double epsilon = 1e-3)
    {
        return fabs(a - b) < epsilon;
    }
};
}

#endif // MATCHER2D_HPP
