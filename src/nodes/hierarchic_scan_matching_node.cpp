#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>

#include <ndt/conversion/convert.hpp>
#include <ndt/data/laserscan.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/visualization/multi_grid.hpp>

struct ScanMatcherNode {
    typedef ndt::matching::MultiGridMatcher2D   MultiGridMatcher2D;
    typedef ndt::visualization::MultiGrid2D     MultiGrid2D;
    typedef ndt::visualization::Point2D         Point2D;
    typedef pcl::PointCloud<pcl::PointXYZ>      PCLPointCloud3D;

    ros::NodeHandle     nh;
    ros::Subscriber     sub;
    ros::Publisher      pub_pcl;
    ros::Publisher      pub_distr;

    ndt::data::LaserScan::Ptr   src;
    MultiGrid2D::ResolutionType resolution_coarse;
    MultiGrid2D::ResolutionType resolution_fine;

    ScanMatcherNode() :
        nh("~"),
        resolution_coarse{2.0, 2.0},
        resolution_fine{1.0, 1.0}
    {
        std::string topic_scan = "/scan";
        std::string topic_pcl  = "/matched";
        std::string topic_distr = "/distribution";
        nh.getParam("topic_scan", topic_scan);
        nh.getParam("topic_pcl", topic_pcl);
        nh.getParam("topic_distr", topic_distr);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic_scan, 1, &ScanMatcherNode::laserscan, this);
        pub_pcl = nh.advertise<PCLPointCloud3D>(topic_pcl, 1);
        pub_distr = nh.advertise<nav_msgs::OccupancyGrid>(topic_distr, 1);

    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        ndt::data::LaserScan dst;
        ndt::conversion::convert(msg, dst);
        if(!src) {
            src.reset(new ndt::data::LaserScan(dst));
        } else {
            MultiGridMatcher2D::Parameters    params_coarse;
            params_coarse.resolution = resolution_coarse;
            MultiGridMatcher2D                matcher_coarse;
            MultiGridMatcher2D::TransformType prior_transform;
            double score_coarse = matcher_coarse.match(dst, *src, prior_transform);

            if(score_coarse < 1e-3)
                prior_transform = MultiGridMatcher2D::TransformType::Identity();

            MultiGridMatcher2D matcher_fine;
            MultiGridMatcher2D::TransformType transform_fine;
            double score_fine = matcher_fine.match(dst, *src, transform_fine, prior_transform);

            PCLPointCloud3D output;

            for(std::size_t i = 0 ; i < src->size ; ++i) {
                Point2D p_bar = transform_fine * src->points[i];
                pcl::PointXYZ pcl_p;
                pcl_p.x = p_bar(0);
                pcl_p.y = p_bar(1);
                output.push_back(pcl_p);
            }

            ndt::data::LaserScan::PointType range = src->range();
            MultiGridMatcher2D::SizeType    size = {std::size_t(range(0) / resolution_fine[0]),
                                                    std::size_t(range(1) / resolution_fine[1])};

            cv::Mat distribution(500,500, CV_8UC3, cv::Scalar());
            MultiGrid2D grid(size, resolution_fine, src->min);
            grid.add(dst);
            ndt::visualization::renderMultiGrid(grid, src->min, src->max, distribution);
            /// output display
            if(score_fine < 0.1){
                std::cout << "-------------------------------" << std::endl;
                std::cout << score_coarse << std::endl;
                std::cout << prior_transform.translation() << std::endl;
                std::cout << prior_transform.rotation() << std::endl;
                std::cout << "-------------------------------" << std::endl;
                std::cout << score_fine << std::endl;
                std::cout << transform_fine.translation() << std::endl;
                std::cout << transform_fine.rotation() << std::endl;
                std::cout << "-------------------------------" << std::endl;
            }
            cv::cvtColor(distribution,  distribution, CV_BGR2GRAY);
            distribution *= 100.0 / 255.0;

            nav_msgs::OccupancyGrid distr_msg;
            distr_msg.header = msg->header;
            distr_msg.info.height = distribution.rows;
            distr_msg.info.width = distribution.cols;
            distr_msg.info.origin.position.x = dst.min(0);
            distr_msg.info.origin.position.y = dst.min(1);
            distr_msg.info.resolution = range(1) / distribution.rows;
            for(int i = 0 ; i < distribution.rows ; ++i) {
                for(int j = 0 ; j < distribution.cols ; ++j)
                    distr_msg.data.push_back(distribution.at<uchar>(distribution.rows - 1 - i, j));
            }

            output.header = pcl_conversions::toPCL(msg->header);
            pub_pcl.publish(output);
            pub_distr.publish(distr_msg);
            src.reset(new ndt::data::LaserScan(dst));
       }
   }
};




int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_hierarchic_grid_matcher_node");
    ScanMatcherNode node;
    ros::spin();

    return 0;
}
