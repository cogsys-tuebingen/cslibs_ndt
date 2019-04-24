#ifndef CSLIBS_NDT_3D_ICP_HPP
#define CSLIBS_NDT_3D_ICP_HPP

#include <cslibs_math_3d/linear/pointcloud.hpp>
#include <cslibs_ndt_3d/matching/icp_params.hpp>
#include <cslibs_ndt_3d/matching/icp_result.hpp>

namespace cslibs_ndt_3d {
namespace matching {
namespace impl {
struct icp {
inline static void apply(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                         const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                         const ParametersWithICP                      &params,
                         const cslibs_math_3d::Transform3d            &initial_transform,
                         ResultWithICP                                &r)
{
    const cslibs_math_3d::Pointcloud3d::points_t &src_points = src->getPoints();
    const cslibs_math_3d::Pointcloud3d::points_t &dst_points = dst->getPoints();
    const std::size_t src_size = src_points.size();
    const std::size_t dst_size = dst_points.size();

    auto sq = [](const double x) {return x * x;};

    const double trans_eps = sq(params.translationEpsilon());
    const double rot_eps = sq(params.rotationEpsilon());
    const double max_distance = sq(params.maxDistanceICP());
    const std::size_t max_iterations = params.maxIterationsICP();

    cslibs_math_3d::Transform3d &transform = r.transform();
    transform = initial_transform;
    cslibs_math_3d::Pointcloud3d::points_t src_points_transformed(src_size);

    std::vector<std::size_t> indices(src_size, std::numeric_limits<std::size_t>::max());
    std::size_t assigned = 0u;

    auto is_assigned = [](const std::size_t index)
    {
        return index < std::numeric_limits<std::size_t>::max();
    };

    cslibs_math_3d::Point3d dst_mean;
    for(const cslibs_math_3d::Point3d &p : dst_points) {
        dst_mean += p;
    }
    dst_mean /= static_cast<double>(dst_size);

    Eigen::Matrix3d &S = r.icpCovariance();

    for(std::size_t i = 0 ; i < max_iterations ; ++i) {
        std::fill(indices.begin(), indices.end(), std::numeric_limits<std::size_t>::max());
        assigned = 0u;

        cslibs_math_3d::Point3d src_mean;

        /// associate
        for(std::size_t s = 0 ; s < src_size ; ++s) {
            cslibs_math_3d::Point3d &sp = src_points_transformed[s];
            sp = transform * src_points[s];
            std::size_t &index = indices[s];
            src_mean += sp;

            double      min_distance = std::numeric_limits<double>::max();
            for(std::size_t d = 0 ; d < dst_size ; ++d) {
                const cslibs_math_3d::Point3d &dp = dst_points[d];
                const double dist = cslibs_math::linear::distance2(dp, sp);
                if(dist < min_distance &&
                        dist < max_distance) {
                    index = d;
                    min_distance = dist;
                }
            }
            assigned += min_distance < max_distance ? 1u : 0u;
        }
        src_mean /= static_cast<double>(src_size);

        S = Eigen::Matrix3d::Zero();
        for(std::size_t s = 0 ; s < src_size ; ++s) {
            const cslibs_math_3d::Point3d &sp = src_points_transformed[s];
            const std::size_t           index = indices[s];
            if(is_assigned(index)) {
                const cslibs_math_3d::Point3d &dp = dst_points[index];
                S += (sp - src_mean).data() * (dp - dst_mean).data().transpose();
            }

        }

        Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3> > svd (S, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d R =(svd.matrixU() * svd.matrixV().transpose()).transpose();
        Eigen::Quaterniond qe(R);

        cslibs_math_3d::Quaterniond   q(qe.x(), qe.y(), qe.z(), qe.w());
        cslibs_math_3d::Transform3d  dt(dst_mean - q * src_mean,
                        q);
        transform *= dt;

        if(dt.translation().length2() < trans_eps ||
                sq(q.angle(cslibs_math_3d::Quaterniond())) < rot_eps) {
            r.icpIterations()  = i;
            r.icpTermination() = ICPTermination::DELTA_EPS;
            return;
        }
        if(static_cast<double>(assigned) / static_cast<double>(src_size)
           < params.minAssignedPoints()) {
            r.iterations() = i;
            r.icpTermination() = ICPTermination::ASSIGNMENT_SUCCESS;
        }

    }

    r.icpIterations() = params.maxIterations();
    r.icpTermination() = ICPTermination::MAX_ITERATIONS;
}
};
}
}
}

#endif // CSLIBS_NDT_3D_ICP_HPP
