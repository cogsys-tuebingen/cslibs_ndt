#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include "../data/laserscan.hpp"
#include "../ndt/matcher2D.hpp"
#include "../convert/convert.hpp"
#include "../tests/optional/visualize.hpp"
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>

struct ScanMatcherNode {
    typedef ndt::NDTMultiGrid<2>    NDTGridType;
    typedef ndt::NDTMatcher2D       NDTMatcher;
    typedef NDTGridType::Resolution NDTResolution;
    typedef pcl::PointCloud<pcl::PointXYZ>
            PCLPointCloudType;

    ros::NodeHandle     nh;
    ros::Subscriber     sub;
    ros::Publisher      pub_pcl;
    ros::Publisher      pub_distr;

    ndt::data::LaserScan::Ptr src;
    NDTResolution             resolution;

    ScanMatcherNode() :
        nh("~"),
        resolution{1.0, 1.0}
    {
        std::string topic_scan = "/scan";
        std::string topic_pcl  = "/matched";
        std::string topic_distr = "/distribution";
        nh.getParam("topic_scan", topic_scan);
        nh.getParam("topic_pcl", topic_pcl);
        nh.getParam("topic_distr", topic_distr);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic_scan, 1, &ScanMatcherNode::laserscan, this);
        pub_pcl = nh.advertise<PCLPointCloudType>(topic_pcl, 1);
        pub_distr = nh.advertise<nav_msgs::OccupancyGrid>(topic_distr, 1);

    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        ndt::data::LaserScan dst;
        ndt::convert::convert(msg, dst, true);
        if(!src) {
            src.reset(new ndt::data::LaserScan(dst));
        } else {
            NDTMatcher matcher(resolution);
            NDTMatcher::TransformType transform;
            double score = matcher.match(*src, dst, transform, 1000, 1e-6, 1e-6);
            transform = NDTMatcher::RotationType(-M_PI_2) * transform;
            if(score < 80) {
                src.reset(new ndt::data::LaserScan(dst));
                return;
            }

            PCLPointCloudType output;

            for(std::size_t i = 0 ; i < dst.size ; ++i) {
                NDTMatcher::PointType p_bar = transform * dst.points[i];
                pcl::PointXYZ pcl_p;
                pcl_p.x = p_bar(0);
                pcl_p.y = p_bar(1);
                output.push_back(pcl_p);
            }

            ndt::data::LaserScan::PointType range = src->range();
            ndt::NDTMultiGrid2D::Size size = {std::size_t(range(0) / resolution[0]),
                                              std::size_t(range(1) / resolution[1])};

            cv::Mat distribution(500,500, CV_8UC3, cv::Scalar());
            ndt::NDTMultiGrid2D grid(size, resolution, src->min);
            grid.add(*src);
            ndt::renderNDTGrid(grid, src->min, src->max, distribution);
            cv::cvtColor(distribution,  distribution, CV_BGR2GRAY);
            distribution *= 0.5;

            nav_msgs::OccupancyGrid distr_msg;
            distr_msg.header = msg->header;
            distr_msg.info.height = distribution.rows;
            distr_msg.info.width = distribution.cols;
            distr_msg.info.origin.position.x = src->min(0);
            distr_msg.info.origin.position.y = -src->min(1);
            distr_msg.info.origin.orientation = tf::createQuaternionMsgFromYaw(-M_PI_2);
            distr_msg.info.resolution = range(1) / distribution.rows;
            for(int i = 0 ; i < distribution.rows * distribution.cols ; ++i)
                distr_msg.data.push_back(distribution.at<uchar>(i));

            output.header = pcl_conversions::toPCL(msg->header);
            pub_pcl.publish(output);
            pub_distr.publish(distr_msg);
            src.reset(new ndt::data::LaserScan(dst));
       }
   }
};




int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_matcher");
    ScanMatcherNode node;
    ros::spin();

    return 0;
}
