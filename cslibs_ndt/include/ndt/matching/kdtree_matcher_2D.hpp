#ifndef KDTREE_MATCHER_2D_HPP
#define KDTREE_MATCHER_2D_HPP

#include <ndt/tree/kdtree.hpp>
#include <ndt/matching/matcher.hpp>
#include <ndt/data/pointcloud.hpp>
#include <ndt/math/angle.hpp>
#include <ndt/math/hausdorff.hpp>

#include <array>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

namespace ndt {
namespace matching {
class KDTreeMatcher2D : public Matcher<2> {



public:

    typedef Eigen::Matrix<double,3,3>               HessianType;
    typedef Eigen::Translation2d                    TranslationType;
    typedef Eigen::Rotation2Dd                      RotationType;
    typedef Eigen::Vector3d                         GradientType;
    typedef PointCloudType::PointType               PointType;
    typedef ndt::math::Distribution<2, true>        DistributionType;
    typedef DistributionType::MatrixType            CovarianceMatrixType;

    typedef ndt::tree::Index<2>                                 IndexType;
    typedef ndt::tree::NodeData<2>                              NodeDataType;
    typedef kdtree::unbuffered::KDTree<IndexType, NodeDataType> KDTreeType;
    typedef KDTreeType::NodeType                                KDNodeType;

    KDTreeMatcher2D(const Parameters &_params = Parameters(),
                    const bool _shuffle = false) :
        BaseClass(_params),
        index(_params.resolution),
        shuffle(_shuffle),
        rotation(0.0)
    {
    }

    virtual ~KDTreeMatcher2D()
    {
    }

    inline typename KDTreeType::Ptr getTree()
    {
        return tree;
    }

    inline double match(const PointCloudType &_dst,
                        const PointCloudType &_src,
                        TransformType        &_transformation,
                        const TransformType  &_prior_transformation = TransformType::Identity()) override
    {
        /// build the ndt grid for the src cloud
        PointCloudType dst = _dst;
        if(shuffle) {
            /// TODO: Move Shuffle to tree insertion!
            std::random_shuffle(dst.points_data.begin(), dst.points_data.end());
        }

        for(auto &p : dst.points_data) {
            tree->insert_bulk(index.create(p), p);
        }
        tree->load_bulk();

        /// reset all the members
        tx        = 0.0;
        ty        = 0.0;
        phi       = 0.0;
        prev_tx   = tx;
        prev_ty   = ty;
        prev_phi  = phi;
        /// gradient and stuff
        /// need 4 fields for that
        /// only solve the maximal score
        max_score             = std::numeric_limits<double>::lowest();
        iteration             = 0;
        lambda                = params.lambda;
        step_corrections      = 0;

        bool converged = false;
        while(!converged) {
            rotation        = RotationType(phi);
            translation     = TranslationType(tx, ty);
            transformation  = translation * rotation * _prior_transformation;

            gradient = GradientType::Zero();
            hessian  = HessianType::Zero();
            score = 0.0;
            double  sin_phi;
            double  cos_phi;
            sincos(phi, &sin_phi, &cos_phi);

            /// calculate the hessian and the gradient for each grid
            for(std::size_t i = 0 ; i < _src.size ; ++i) {
                if(_src.mask[i] == PointCloudType::VALID) {
                    PointType p = transformation * _src.points[i];

                    KDNodeType *node = tree->find(index.create(p));
                    if(node == nullptr)
                        continue;
                    if(!node->data.distribution)
                        continue;

                    DistributionType &distribution = *(node->data.distribution);
                    if(distribution.getN() < 3)
                        continue;

                    PointType            mean;
                    PointType            q;
                    CovarianceMatrixType inverse_covariance;
                    PointType            jac = PointType::Zero();
                    PointType            hes = PointType::Zero();

                    double s = distribution.sampleNonNormalized(p, q);
                    distribution.getMean(mean);
                    distribution.getInverseCovariance(inverse_covariance);

                    /// at this point, s must be greater than 0.0, since we work with non-normalized Gaussians.
                    /// if s is 0.0 we do no need to count the sample in
                    if(s > 1e-3) {
                        score += s;
                        const PointType q_inverse_covariance = q.transpose() * inverse_covariance;

                        jac(0) = -p(0) * sin_phi - p(1) * cos_phi;
                        jac(1) =  p(0) * cos_phi - p(1) * sin_phi;
                        hes(0) = -p(0) * cos_phi + p(1) * sin_phi;
                        hes(1) = -p(0) * sin_phi - p(1) * cos_phi;

                        //                            /// gradient computation
                        double g_dot = q_inverse_covariance.transpose() * jac;
                        gradient(0) += s * q_inverse_covariance(0);
                        gradient(1) += s * q_inverse_covariance(1);
                        gradient(2) += s * g_dot;

                        hessian(0,0)+= -s * (
                                    (-q_inverse_covariance(0) * -q_inverse_covariance(0))   /// (1)
                                    +(-inverse_covariance(0,0)));                           /// (3)
                        hessian(1,0)+= -s * (
                                    (-q_inverse_covariance(1) * -q_inverse_covariance(0))   /// (1)
                                    +(-inverse_covariance(1,0)));                           /// (3)
                        hessian(2,0)+= -s * (
                                    (-g_dot * -q_inverse_covariance(0))                     /// (1)
                                    +(-inverse_covariance.row(0).dot(jac)));                /// (3)
                        hessian(0,1)+= -s * (
                                    (-q_inverse_covariance(0) * -q_inverse_covariance(1))   /// (1)
                                    +(-inverse_covariance(0,1)));                           /// (3)
                        hessian(1,1)+= -s * (
                                    (-q_inverse_covariance(1) * -q_inverse_covariance(1))   /// (1)
                                    +(-inverse_covariance(1,1)));                           /// (3)
                        hessian(2,1)+= -s * (
                                    (-g_dot * -q_inverse_covariance(1))                     /// (1)
                                    +(-inverse_covariance.row(1).dot(jac)));                /// (3)
                        hessian(0,2)+= -s * (
                                    (-q_inverse_covariance(0) * -g_dot)                     /// (1)
                                    +(-jac.transpose() * inverse_covariance.col(0)));        /// (3)
                        hessian(1,2)+= -s * (
                                    (-q_inverse_covariance(1) * -g_dot )                    /// (1)
                                    +(-jac.transpose() * inverse_covariance.col(1)));        /// (3)
                        hessian(2,2)+= -s * (
                                    (-g_dot * -g_dot)                                       /// (1)
                                    +(-q_inverse_covariance.transpose() * hes)              /// (2)
                                    +(-jac.transpose() * inverse_covariance * jac));        /// (3)


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

            /// now we have to check wether the score increased or not
            /// if not, we have to adjust the step size
            if(score < max_score) {
                tx      = prev_tx;
                ty      = prev_ty;
                phi     = prev_phi;
                lambda *= params.alpha;
                rotation        = RotationType(phi);
                translation     = TranslationType(tx, ty);
                transformation  = translation * rotation * _prior_transformation;
                ++step_corrections;
            } else {
                if(iteration > 0 &&
                        epsTrans(tx, prev_tx) &&
                        epsTrans(ty, prev_ty) &&
                        epsRot(phi, prev_phi)) {
                    break;
                }

                max_score = score;
                prev_tx   = tx;
                prev_ty   = ty;
                prev_phi  = phi;
                step_corrections = 0;
//                lambda = params.lambda;
                if(score > max_score)
                    lambda /= params.alpha;
            }

            if(step_corrections >= params.max_step_corrections) {
                break;
            }

            if(iteration >= params.max_iterations) {
                break;
            }


            /// positive definiteness
            Eigen::EigenSolver<HessianType> hessian_solver(hessian, false);
            double min_eigen_value = hessian_solver.eigenvalues().real().minCoeff();
            if(min_eigen_value < 0) {
                double lambda = 1.1 * min_eigen_value - 1;
                hessian += Eigen::Vector3d (-lambda, -lambda, -lambda).asDiagonal();
            }

            /// solve equeation here
            delta_p = GradientType::Zero();
            delta_p = (-hessian).fullPivLu().solve(gradient);
            tx  += delta_p(0) * lambda(0);
            ty  += delta_p(1) * lambda(1);
            phi += delta_p(2) * lambda(2);
            phi = math::wrapAngle(phi);

            ++iteration;
        }

        rotation        = RotationType(phi);
        translation     = TranslationType(tx, ty);
        transformation  = translation * rotation * _prior_transformation;
        _transformation = transformation;

        /// todo exchange through kdtree nn ...
        std::size_t accepted = 0;
        std::size_t size_valid = 0;
        for(std::size_t i = 0 ; i < _dst.size ; ++i) {
            if(_src.mask[i] == PointCloudType::VALID) {
                PointType p = transformation * _src.points[i];
                double h = ndt::math::hausdorff<2>(p, _dst);
                if(h < 0.1)
                    ++accepted;
                ++size_valid;
            }
        }

        return accepted / (double) size_valid;

    }

    void printDebugInfo()
    {
        std::cout << "----------------------------------" << std::endl;
        std::cout << "score : " << max_score << std::endl;
        std::cout << "translation : " << std::endl;
        std::cout << transformation.translation() << std::endl;
        std::cout << "rotation : " << std::endl;
        std::cout << transformation.rotation() << std::endl;
        std::cout << "iterations : " << iteration << std::endl;
        std::cout << "lambda : " << lambda << std::endl;
        std::cout << "step corrections : " << step_corrections << std::endl;
        std::cout << "----------------------------------" << std::endl;

    }

private:
    KDTreeType::Ptr             tree;
    IndexType                   index;
    bool                        shuffle;

    double tx;
    double ty;
    double phi;
    double prev_tx;
    double prev_ty;
    double prev_phi;
    double max_score;

    TransformType    transformation;
    TranslationType  translation;
    RotationType     rotation;

    double           score;
    GradientType     gradient;
    HessianType      hessian;
    GradientType     delta_p;

    std::size_t      iteration;
    LambdaType       lambda;
    double           alpha;
    std::size_t      step_corrections;

};
}
}
#endif // KDTREE_MATCHER_2D_HPP
