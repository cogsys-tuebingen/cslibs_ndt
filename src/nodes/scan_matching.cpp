#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include "../data/laserscan.hpp"
#include "../ndt/ndt_multi_grid.h"
#include <eigen3/Eigen/Jacobi>

namespace ndt {
struct ScanVisualizerNode {
    ros::NodeHandle   nh;
    ros::Subscriber   sub;

    double            resolution;
    double            size;
    double            margin;

    std::size_t       conv_iter;
    double            conv_eps;

    NDTMultiGrid::Ptr ndt_grid;

    ScanVisualizerNode() :
        nh("~"),
        resolution(1.0),
        size(0.0),
        margin(0.1),
        conv_iter(100),
        conv_eps(1e-3)
    {
        std::string topic("/scan");
        nh.getParam("topic", topic);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic, 1, &ScanVisualizerNode::laserscan, this);

    }

    void laserscan(const sensor_msgs::LaserScan::ConstPtr &msg)
    {
        LaserScan scan(msg);

        if(!ndt_grid) {
            size = msg->range_max;
            std::cout << size << std::endl;
            ndt_grid.reset(new NDTMultiGrid(size + margin,
                                            size + margin,
                                            resolution));
            insert(scan, ndt_grid);
        } else {
            for(int i = 0 ; i < 3 ; ++i) {
                NDTGrid &g = ndt_grid->at(i);
                for(std::size_t j = 0 ; j < g.dim_x() ; ++j) {
                    for(std::size_t k = 0 ; k < g.dim_y() ; ++k) {
                        RollingDistribution &r = g.at(j,k);
                        if(r.n() >= 3) {
                            Point      mean;
                            Covariance cov;
                            r.mean(mean);
                            r.covariance(cov);
                            std::cout << "[ " << i << " " << j << " " << k  << " ]" << std::endl;
                            std::cout << mean << std::endl;
                            std::cout << cov << std::endl;
                        }
                    }
                }

            }
//            if(size != msg->range_max) {
//                std::cerr << "Laserscan properties are not matching!" << std::endl;
//                return;
//            }

//            TransformMatrix m;
//            match(scan, m);

//            std::cout << m << std::endl;

//            insert(scan, ndt_grid);
        }
    }

    void insert(const LaserScan   &scan,
                NDTMultiGrid::Ptr &ndt_grid)
    {
        for(std::size_t i = 0 ; i < scan.size ; ++i) {
            if(scan.mask[i] == LaserScan::VALID) {
                ndt_grid->add(scan.points[i]);
            }
        }
    }

    void match(const LaserScan &scan,
               TransformMatrix &transform)
    {
        /// 1. Build the Normal Distribution Transform of the first scan
        ndt_grid.reset(new NDTMultiGrid(size + margin,
                                        size + margin,
                                        resolution));
        insert(scan, ndt_grid);

//        /// 2. Intialize the estimate for the paramters
//        transform = TransformMatrix::Identity();

//        double theta = 0.0;
//        double tx = 0.0;
//        double ty = 0.0;

//        Point *points = new Point[scan.size];
//        double score = 0.0;
//        double          s = 0.0;
//        PointTransposed q_cov = Point::Zero();
//        Transform       g     = Transform::Zero();
//        TransformMatrix mp    = TransformMatrix::Identity();
//        Hessian         h     = Hessian::Zero();
//        Point           mean;
//        Covariance      inv_cov;

//        Jacobian jac = Jacobian::Identity();
//        Point    hes = Point::Zero();

//        bool converged = false;
//        std::size_t iteration = 0;
//        while(!converged && iteration < conv_iter) {
//            Eigen::Matrix2d r;
//            r(0,0) =  cos(theta);
//            r(0,1) = -sin(theta);
//            r(1,0) =  sin(theta);
//            r(1,1) =  cos(theta);
//            Eigen::Vector2d t;
//            t(0) = tx;
//            t(1) = ty;

//            /// 3. For each sample fo the second scan: Map the reconstructed 2D point into the
//            ///    the coordinate frame of the first scan according to the parameters.
//            /// 4. Determine the corresponding normal distributions for each mapped point.
//            /// 5. The score for the parameters is determined by evaluation the distribution
//            ///    for each point and summing the result.
//            /// 6. Calculate the new parameter estimate by traying to optimize the score.
//            ///    This is done by performin one step of Newton's Algorithm.
//            /// 7. Go to 3. until a convergence criterion is met.

//            score = 0.0;
//            for(std::size_t i = 0; i < scan.size; ++i) {
//                if(scan.mask[i] == LaserScan::VALID) {
//                    points[i] = r * scan.points[i] + t;
//                    double x = points[i](0); double y = points[i](1);
//                    double sin, cos;
//                    sincos(theta, &sin, &cos);
//                    jac(0,2) = -sin * x -cos * y;
//                    jac(1,2) =  cos * x -sin * y;
//                    hes(0)   = -cos * x +sin * y;
//                    hes(1)   = -sin * x -cos * y;

//                    for(int j = 0 ; j < 4 ; ++j) {
//                        ndt_grid->at(j).sample(points[i], s, mean, inv_cov);
//                        q_cov = (points[i] - mean).transpose() * inv_cov;
//                        for(std::size_t k = 0 ; k < 3 ; ++k)
//                            g(k)   += double(q_cov * jac.col(k)) * s;
//                        for(std::size_t k = 0 ; k < 3 ; ++k) {
//                            for(std::size_t l = 0 ; l < 3; ++l) {
//                                h(k,l) += double((q_cov  * jac.col(k)))
//                                        * double((-q_cov * jac.col(l)))
//                                        + double((jac.col(l).transpose() * inv_cov * jac.col(k)));
//                                h(k,l) *= s;
//                            }
//                        }
//                        h(2,2) += double(q_cov * hes) * s;
//                        score += s;
//                    }
//                }
//            }

            /// check for positive definiteness
//            Eigen::EigenSolver<Hessian> solver;
//            solver.compute(h);
//            Eigen::Vector3cd eigen_values_complex  = solver.eigenvalues();
//            if(eigen_values_complex(0).real() < 0.0 ||
//                    eigen_values_complex(1).real() < 0.0 ||
//                        eigen_values_complex(2).real() < 0.0) {
//                double max = std::numeric_limits<double>::min();
//                for(std::size_t i = 0 ; i < 3; ++i) {
//                    for(std::size_t j = 0 ; j < 3 ; ++j) {
//                        double a = fabs(h(i,j));
//                        if(a > max)
//                            max = a;
//                    }
//                }
//                h += Eigen::Matrix3d::Identity() * (1.5 * max);
//            }

            /// solve the equation
            /// save the complete transformation
//            Transform       p;
//            Eigen::JacobiSVD<Hessian> svd(h, Eigen::ComputeFullU | Eigen::ComputeFullV);
//            p = svd.solve(-g);
//            mp(0,0) = cos(p(2));
//            mp(0,1) = -sin(p(2));
//            mp(1,0) = sin(p(2));
//            mp(1,1) = cos(p(2));
//            mp(0,2) = p(0);
//            mp(1,2) = p(1);
//            transform = mp * transform;

//            if(fabs(theta) < conv_eps &&
//                    fabs(tx) < conv_eps &&
//                    fabs(ty) < conv_eps) {
//                break;
//            }

//            ++iteration;
//        }
//        delete[] points;
     }




};
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "scan_matcher");
    ndt::ScanVisualizerNode sn;

    ros::spin();

    return 0;
}
