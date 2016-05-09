#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>
#include <chrono>

#include <ndt/conversion/convert_ros.hpp>
#include <ndt/data/laserscan.hpp>
#include <ndt/matching/multi_grid_matcher_2D.hpp>
#include <ndt/visualization/multi_grid.hpp>

struct ScanMatcherNode {
    typedef ndt::matching::MultiGridMatcher2D   MatcherType;
    typedef ndt::visualization::MultiGrid2D     MultiGrid2D;
    typedef ndt::visualization::Point2D         Point2D;
    typedef pcl::PointCloud<pcl::PointXYZ>      PCLPointCloudType;

    ros::NodeHandle       nh;
    ros::Subscriber       sub;
    ros::Publisher        pub_pcl;
    ros::Publisher        pub_distr;
    tf::TransformListener tf;

    float                       range_min;
    float                       range_max;
    ndt::data::LaserScan::Ptr   src;
    tf::StampedTransform        src_transform;
    MultiGrid2D::ResolutionType resolution;
    std::size_t                 failed;
    std::size_t                 all;
    std::size_t                 count;

    ScanMatcherNode() :
        nh("~"),
        range_min(-1.f),
        range_max(-1.f),
        resolution{1.0, 1.0},
        failed(0),
        all(0),
        count(0)
    {
        std::string topic_scan  = "/scan";
        std::string topic_pcl   = "/matched";
        std::string topic_distr = "/distribution";

        nh.getParam("topic_scan", topic_scan);
        nh.getParam("topic_pcl", topic_pcl);
        nh.getParam("topic_distr", topic_distr);
        nh.getParam("range_min", range_min);
        nh.getParam("range_max", range_max);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic_scan, 1, &ScanMatcherNode::laserscan, this);
        pub_pcl = nh.advertise<PCLPointCloudType>(topic_pcl, 1);
        pub_distr = nh.advertise<nav_msgs::OccupancyGrid>(topic_distr, 1);

        std::cout << "topic_scan : " << topic_scan << std::endl;
        std::cout << "topic_pcl : " << topic_pcl << std::endl;
        std::cout << "topic_distr : " << topic_distr << std::endl;
        std::cout << "range_min : " << range_min << std::endl;
        std::cout << "range_max : " << range_max << std::endl;
    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();

        /// match old points to the current ones
        ndt::data::LaserScan dst;
        ndt::conversion::convert(msg, dst, range_min, range_max);

        tf::StampedTransform dst_transform;
        try{
            tf.waitForTransform(msg->header.frame_id, "/odom", msg->header.stamp, ros::Duration(0.1));
            tf.lookupTransform(msg->header.frame_id, "/odom", msg->header.stamp, dst_transform);
        }
        catch (const tf::TransformException &ex){
            ROS_ERROR("%s",ex.what());
            return;
        }

        if(!src) {
            src_transform = dst_transform;
            src.reset(new ndt::data::LaserScan(dst));
        } else {
            tf::Transform  diff = (dst_transform.inverse() * src_transform);

            MatcherType::Parameters params;
            params.max_step_corrections = 5;
            params.eps_rot = 1e-6;
            params.lambda = MatcherType::LambdaType::Constant(1.0);
            MatcherType matcher(params);
            MatcherType::RotationType    rotation(tf::getYaw(diff.getRotation()));
            MatcherType::TranslationType translation(diff.getOrigin().x(), diff.getOrigin().y());
            MatcherType::TransformType   prior_transform = translation * rotation;
            MatcherType::TransformType   transform = prior_transform;
            double score = matcher.match(dst, *src, transform);

            std::chrono::microseconds elapsed =
                    std::chrono::duration_cast<std::chrono::microseconds>
                    (std::chrono::system_clock::now() - start);
            std::cout << "elapsed " << elapsed.count() / 1000.0 << " ms" << std::endl;

            PCLPointCloudType output;
            for(std::size_t i = 0 ; i < src->size ; ++i) {
                Point2D p_bar = transform * src->points[i];
                pcl::PointXYZ pcl_p;
                pcl_p.x = p_bar(0);
                pcl_p.y = p_bar(1);
                output.push_back(pcl_p);
            }

//            ndt::data::LaserScan::PointType range = src->range();
//            MultiGrid2D::SizeType  size = {std::size_t(range(0) / resolution[0]),
//                                           std::size_t(range(1) / resolution[1])};

//            cv::Mat distribution(500,500, CV_8UC3, cv::Scalar());
//            MultiGrid2D grid(size, resolution, dst.min);
//            grid.add(dst);
//            ndt::visualization::renderMultiGrid(grid,
//                                                dst.min,
//                                                dst.max,
//                                                distribution);
            matcher.printDebugInfo();
            /// output display
            if(score < 50){
//                std::cout << "-------------------------------" << std::endl;
//                std::cout << score << std::endl;
//                std::cout << transform.translation() << std::endl;
//                std::cout << transform.rotation() << std::endl;
//                std::cout << "-------------------------------" << std::endl;

//                {
//                    std::stringstream ss_dst;
//                    std::stringstream ss_src;
//                    ss_dst << "/tmp/dst_" << failed << ".scan";
//                    ss_src << "/tmp/src_" << failed << ".scan";
//                    dst.save(ss_dst.str());
//                    src->save(ss_src.str());
//                }
                ++failed;
            }
//            cv::cvtColor(distribution,  distribution, CV_BGR2GRAY);
//            distribution *= 100.0 / 255.0;

//            nav_msgs::OccupancyGrid distr_msg;
//            distr_msg.header = msg->header;
//            distr_msg.info.height = distribution.rows;
//            distr_msg.info.width = distribution.cols;
//            distr_msg.info.origin.position.x = dst.min(0);
//            distr_msg.info.origin.position.y = dst.min(1);
//            distr_msg.info.resolution = range(1) / distribution.rows;
//            for(int i = 0 ; i < distribution.rows ; ++i) {
//                for(int j = 0 ; j < distribution.cols ; ++j)
//                    distr_msg.data.push_back(distribution.at<uchar>(distribution.rows - 1 - i, j));
//            }

//            pub_distr.publish(distr_msg);

            output.header = pcl_conversions::toPCL(msg->header);
            pub_pcl.publish(output);

            if(count >= 5) {
                src.reset(new ndt::data::LaserScan(dst));
                src_transform = dst_transform;
                count = 0;
            }
            ++count;

            ++all;


            /// std::cout << "success : " << failed << " " << all << " => " << (1.0 - failed / (double) all) * 100.0 << std::endl;
        }
    }
};




int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_multi_grid_matcher_node");
    ScanMatcherNode node;
    ros::spin();

    return 0;
}
