#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>
#include <chrono>

#include <ndt/conversion/convert.hpp>
#include <ndt/data/laserscan.hpp>
#include <ndt/matching/kdtree_matcher_2D.hpp>
#include <ndt/tree/kdtree.hpp>
#include <ndt/visualization/kdtree.hpp>

struct ScanMatcherNode {
    typedef ndt::tree::KDTreeNode<2>        KDTreeNodeType;
    typedef KDTreeNodeType::KDTreeType      KDTreeType;
    typedef ndt::tree::KDTreeInterface<2>   KDTreeInterfaceType;
    typedef KDTreeInterfaceType::ResolutionType ResolutionType;
    typedef pcl::PointCloud<pcl::PointXYZ>  PCLPointCloudType;
    typedef ndt::matching::KDTreeMatcher2D  MatcherType;

    ros::NodeHandle     nh;
    ros::Subscriber     sub;
    ros::Publisher      pub_pcl;
    ros::Publisher      pub_distr;
    tf::TransformListener tf;

    float                     range_min;
    float                     range_max;
    ndt::data::LaserScan::Ptr src;
    tf::StampedTransform      src_transform;
    ResolutionType            resolution;
    std::size_t               failed;
    std::size_t               all;

    ScanMatcherNode() :
        nh("~"),
        range_min(-1.f),
        range_max(-1.f),
        resolution{1.0, 1.0},
        failed(0),
        all(0)
    {
        std::string topic_scan = "/scan";
        std::string topic_pcl  = "/matched";
        std::string topic_distr = "/distribution";
        nh.getParam("topic_scan", topic_scan);
        nh.getParam("topic_pcl", topic_pcl);
        nh.getParam("topic_distr", topic_distr);
        nh.getParam("range_min", range_min);
        nh.getParam("range_max", range_max);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic_scan, 1, &ScanMatcherNode::laserscan, this);
        pub_pcl = nh.advertise<PCLPointCloudType>(topic_pcl, 1);
        pub_distr = nh.advertise<nav_msgs::OccupancyGrid>(topic_distr, 1);

    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();

        ndt::data::LaserScan dst;
        ndt::conversion::convert(msg, dst);

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
            src.reset(new ndt::data::LaserScan(dst));
        } else {
            tf::Transform  diff = dst_transform.inverse() * src_transform;

            MatcherType matcher;
            MatcherType::RotationType    rotation(tf::getYaw(diff.getRotation()));
            MatcherType::TranslationType translation(diff.getOrigin().x(), diff.getOrigin().y());
            MatcherType::TransformType   transform = translation * rotation;
            double score = matcher.match(dst, *src, transform);


            std::chrono::duration<double> elapsed =
                    std::chrono::system_clock::now() - start;
            std::cout << "elapsed " << elapsed.count() * 1000.0 << "ms" << std::endl;

            PCLPointCloudType output;

            for(std::size_t i = 0 ; i < src->size ; ++i) {
                MatcherType::PointType p_bar = transform * src->points[i];
                pcl::PointXYZ pcl_p;
                pcl_p.x = p_bar(0);
                pcl_p.y = p_bar(1);
                output.push_back(pcl_p);
            }

            if(score < 0.1){
                std::cout << "-------------------------------" << std::endl;
                std::cout << score << std::endl;
                std::cout << transform.translation() << std::endl;
                std::cout << transform.rotation() << std::endl;
                std::cout << "-------------------------------" << std::endl;
                ++failed;
            }

            KDTreeType::Ptr     tree;
            KDTreeInterfaceType tree_interface(resolution);
            tree_interface.insert(dst, tree);

            cv::Mat distribution(500,500, CV_8UC3, cv::Scalar());
            ndt::visualization::renderTree(tree, tree_interface, dst.min, dst.max, distribution);

            cv::cvtColor(distribution,  distribution, CV_BGR2GRAY);
            distribution *= 100.0 / 255.0;

            nav_msgs::OccupancyGrid distr_msg;
            distr_msg.header = msg->header;
            distr_msg.info.height = distribution.rows;
            distr_msg.info.width = distribution.cols;
            distr_msg.info.origin.position.x = dst.min(0);
            distr_msg.info.origin.position.y = dst.min(1);
            distr_msg.info.resolution = dst.range()(1) / distribution.rows;
            for(int i = 0 ; i < distribution.rows ; ++i) {
                for(int j = 0 ; j < distribution.cols ; ++j)
                    distr_msg.data.push_back(distribution.at<uchar>(distribution.rows - 1 - i, j));
            }

            output.header = pcl_conversions::toPCL(msg->header);
            pub_pcl.publish(output);
            pub_distr.publish(distr_msg);
            src.reset(new ndt::data::LaserScan(dst));
            ++all;

            src_transform = dst_transform;
            std::cout << "success : " << failed << " " << all << " => " << (1.0 - failed / (double) all) * 100.0 << std::endl;
       }
   }
};




int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_kdtree_matcher_node");
    ScanMatcherNode node;
    ros::spin();

    return 0;
}
