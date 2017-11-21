#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl_ros/point_cloud.h>
#include <nav_msgs/OccupancyGrid.h>
#include <tf/tf.h>
#include <fstream>
#include <chrono>

#include <ndt/conversion/convert_ros.hpp>
#include <ndt/data/laserscan.hpp>
#include <ndt/matching/kdtree_matcher_2D.hpp>
#include <ndt/tree/kdtree.hpp>
#include <ndt/visualization/kdtree.hpp>
#include <kdtree/kdtree_clustering.hpp>
#include <kdtree/kdtree_dotty.hpp>

struct ScanClusterNode {
    typedef ndt::tree::KDTreeNode<2>            KDTreeNodeType;
    typedef KDTreeNodeType::KDTreeType          KDTreeType;
    typedef ndt::tree::KDTreeInterface<2>       KDTreeInterfaceType;
    typedef KDTreeInterfaceType::ResolutionType ResolutionType;
    typedef pcl::PointCloud<pcl::PointXYZRGB>   PCLPointCloudType;
    typedef ndt::matching::KDTreeMatcher2D      MatcherType;

    ros::NodeHandle     nh;
    ros::Subscriber     sub;
    ros::Publisher      pub_pcl;

    ndt::data::LaserScan::Ptr dst;
    ResolutionType            resolution;

    ScanClusterNode() :
        nh("~"),
        resolution{0.25, 0.25}
    {
        std::string topic_scan = "/scan";
        std::string topic_pcl  = "/matched";
        nh.getParam("topic_scan", topic_scan);
        nh.getParam("topic_pcl", topic_pcl);

        sub = nh.subscribe<sensor_msgs::LaserScan>(topic_scan, 1, &ScanClusterNode::laserscan, this);
        pub_pcl = nh.advertise<PCLPointCloudType>(topic_pcl, 1);

    }

    void laserscan(const sensor_msgs::LaserScanConstPtr &msg)
    {
        std::chrono::time_point<std::chrono::system_clock> start =
                std::chrono::system_clock::now();
        ndt::data::LaserScan scan;
        ndt::conversion::convert(msg, scan);
        KDTreeType::Ptr     tree;
        KDTreeInterfaceType tree_interface(resolution);
                            tree_interface.insert(scan, tree);
                            tree_interface.cluster(tree);
        KDTreeInterfaceType::DistributionMapType distributions;
        tree_interface.getClusterDistributions(tree, distributions);

        std::chrono::duration<double> elapsed =
                std::chrono::system_clock::now() - start;
        std::cout << "elapsed " << elapsed.count() << "s" << std::endl;
        std::cout << "elapsed " << elapsed.count() * 1000.0 << "ms" << std::endl;

        PCLPointCloudType output;
        kdtree::ColorMap<true> color_palette;
        std::vector<KDTreeNodeType::Ptr> leaves;
        tree->getLeaves(leaves);
        for(KDTreeNodeType::Ptr &leaf : leaves) {
            KDTreeNodeType *node = (KDTreeNodeType*) leaf.get();
            kdtree::Color color;
            color_palette.getColor(node->cluster, color);
            for(auto &p : node->points) {
                pcl::PointXYZRGB pcl_p;
                pcl_p.x = p(0);
                pcl_p.y = p(1);
                pcl_p.r = color.value[0] * 255;
                pcl_p.g = (color.value[0] + color.value[1]) * 127;
                pcl_p.b = (color.value[1] + color.value[2]) * 127;
                output.push_back(pcl_p);
            }
        }
        output.header = pcl_conversions::toPCL(msg->header);
        pub_pcl.publish(output);
   }
};




int main(int argc, char *argv[])
{
    ros::init(argc, argv, "ndt_kdtree_cluster_node");
    ScanClusterNode node;
    ros::spin();

    return 0;
}
