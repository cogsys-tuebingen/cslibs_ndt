#ifndef CSLIBS_NDT_3D_ICP_HPP
#define CSLIBS_NDT_3D_ICP_HPP

#include <cslibs_math_3d/linear/pointcloud.hpp>
#include <cslibs_ndt_3d/matching/params.hpp>
#include <cslibs_ndt_3d/matching/result.hpp>

namespace cslibs_ndt_3d {
namespace matching {
namespace impl {
inline void apply(const cslibs_math_3d::Pointcloud3d::ConstPtr &src,
                  const cslibs_math_3d::Pointcloud3d::ConstPtr &dst,
                  const ParametersWithICP                      &params,
                  ResultWithICP                                &r)
{
    const Pointcloud3d::points_t &src_points = src->getPoints();
    const Pointcloud3d::points_t &dst_points = dst->getPoints();
    const std::size_t src_size = src_points.size();
    const std::size_t dst_size = dst_points.size();

    auto sq = [](const double x) {return x * x;};

    const double trans_eps = sq(params.transEps());
    const double rot_eps = sq(params.rotEps());
    const double max_distance = sq(params.maxDistanceICP());
    const std::size_t max_iterations = params.maxIterationsICP();

    Transform3d &transform = r.transform();
    transform = params.transform();
    Pointcloud3d::points_t src_points_transformed(src_size);

//    std::vector<std::size_t> indices(src_size, std::numeric_limits<std::size_t>::max());

//    auto is_assigned = [](const std::size_t index)
//    {
//        return index < std::numeric_limits<std::size_t>::max();
//    };

//    Point3d dst_mean;
//    for(const Point3d &p : dst_points) {
//        dst_mean += p;
//    }
//    dst_mean /= static_cast<double>(dst_size);


//    Eigen::Matrix3d &S = r.covariance();

//    for(std::size_t i = 0 ; i < max_iterations ; ++i) {
//        std::fill(indices.begin(), indices.end(), std::numeric_limits<std::size_t>::max());

//        Point3d src_mean;

//        /// associate
//        for(std::size_t s = 0 ; s < src_size ; ++s) {
//            Point3d &sp = src_points_transformed[s];
//            sp = transform * src_points[s];
//            std::size_t &index = indices[s];
//            src_mean += sp;

//            double      min_distance = std::numeric_limits<double>::max();
//            for(std::size_t d = 0 ; d < dst_size ; ++d) {
//                const Point3d &dp = dst_points[d];
//                const double dist = distance2(dp, sp);
//                if(dist < min_distance &&
//                        dist < max_distance) {
//                    index = d;
//                    min_distance = dist;
//                }
//            }
//        }
//        src_mean /= static_cast<double>(src_size);

//        for(std::size_t s = 0 ; s < src_size ; ++s) {
//            const Point3d &sp = src_points_transformed[s];
//            const std::size_t index = indices[s];
//            if(is_assigned(index)) {
//                const Point3d &dp = dst_points[index];
//                S += (sp - src_mean).data() * (dp - dst_mean).data().transpose();
//            }

//        }

//        Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3> > svd (S, Eigen::ComputeFullU | Eigen::ComputeFullV);
//        Eigen::Matrix3d R =(svd.matrixU() * svd.matrixV().transpose()).transpose();
//        // Eigen::Matrix<double, 3, 1> T = dst_mean.data() - R * src_mean.data();
//        Eigen::Quaterniond qe(R);

//        Quaternion   q(qe.x(), qe.y(), qe.z(), qe.w());
//        Transform3d  dt(dst_mean - q * src_mean,
//                        q);
//        transform *= dt;

//        if(dt.translation().length2() < trans_eps ||
//                sq(q.angle(Quaternion())) < rot_eps) {
//            r.iterations()  = i;
//            r.termination() = Result::EPS;
//            return;
//        }
//    }

//    r.iterations() = params.maxIterations();
//    r.termination() = Result::ITERATIONS;
}
}
}
}

#endif // CSLIBS_NDT_3D_ICP_HPP
