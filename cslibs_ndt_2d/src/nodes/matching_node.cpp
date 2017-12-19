#include <ros/ros.h>
#include <sensor_msgs/LaserScan.h>

#include <tf/tf.h>
#include <tf/transform_listener.h>
#include <tf/transform_datatypes.h>

#include <cslibs_ndt_2d/dynamic_maps/algorithms/matcher.hpp>
#include <cslibs_time/time.hpp>
#include <cslibs_math_ros/sensor_msgs/conversion_2d.hpp>
#include <cslibs_math_2d/linear/transform.hpp>

namespace cslibs_ndt_2d {
class MatchingNode {
public:
    using matcher_t = dynamic_maps::algorithms::Matcher;

    MatchingNode() :
        nh_("~")
    {
        std::string topic_sub, topic_pub_src, topic_pub_match, topic_pub_dst, topic_pub_map;

        nh_.param<std::string>("topic_scan",  topic_sub,       "/scan");
        nh_.param<std::string>("topic_src",   topic_pub_src,   "/ndt_matcher_2d/src");
        nh_.param<std::string>("topic_match", topic_pub_match, "/ndt_matcher_2d/matched");
        nh_.param<std::string>("topic_dst",   topic_pub_dst,   "/ndt_matcher_2d/dst");
        nh_.param<std::string>("topic_map",   topic_pub_map,   "/ndt_matcher_2d/map");
        nh_.param<std::string>("odom_frame",  odom_frame_,     "/odom");

        sub_       = nh_.subscribe<sensor_msgs::LaserScan> (topic_sub,       1, &MatchingNode::laserscan, this);
        pub_src_   = nh_.advertise<sensor_msgs::PointCloud>(topic_pub_src,   1);
        pub_match_ = nh_.advertise<sensor_msgs::PointCloud>(topic_pub_match, 1);
        pub_dst_   = nh_.advertise<sensor_msgs::PointCloud>(topic_pub_dst,   1);
        pub_map_   = nh_.advertise<nav_msgs::OccupancyGrid>(topic_pub_map,   1);
    }

private:
    void laserscan(
            const sensor_msgs::LaserScanConstPtr & msg)
    {
        cslibs_time::Time start = cslibs_time::Time::now();

        // match old point to the current ones -> current laserscan is dst
        cslibs_math_2d::Pointcloud2d::Ptr dst;
        cslibs_math_ros::sensor_msgs::conversion_2d::from(msg, dst);

        // estimate current pose with odom
        tf::StampedTransform dst_transform;
        try {
            tf_.waitForTransform(msg->header.frame_id, odom_frame_, msg->header.stamp, ros::Duration(0.1));
            tf_.lookupTransform(msg->header.frame_id,  odom_frame_, msg->header.stamp, dst_transform);
        } catch (const tf::TransformException & ex) {
            ROS_ERROR("%s", ex.what());
            return;
        }

        // if a previous laserscan was captured, match it to this one
        if (src_) {

            // initial transform estimate, based on odom
            tf::Transform diff = src_transform_ * dst_transform.inverse();
            cslibs_math_2d::Transform2d transform(diff.getOrigin().x(),
                                                  diff.getOrigin().y(),
                                                  tf::getYaw(diff.getRotation()));
            cslibs_math_2d::Transform2d initial_transform(transform);
            std::cout << transform << std::endl;

            // initialize matching parameters and matcher
            matcher_t::Parameters params;
            int max_iterations, max_step_corrections;
            double lambda = 0.1;
            bool use_odom = false;
            nh_.param<double>("resolution",           params.resolution_,   1.0);
            nh_.param<double>("eps_rot",              params.eps_rot_,      1e-3);
            nh_.param<double>("eps_trans",            params.eps_trans_,    1e-3);
            nh_.param<int>   ("max_iterations",       max_iterations,       100);
            nh_.param<int>   ("max_step_corrections", max_step_corrections, 10);
            nh_.param<double>("alpha",                params.alpha_,        2.0);
            nh_.param<double>("lambda",               lambda,               0.1);
            nh_.param<bool>  ("use_odom",             use_odom,             false);
            params.max_iterations_       = max_iterations;
            params.max_step_corrections_ = max_step_corrections;
            params.lambda_               = matcher_t::lambda_t::Constant(lambda);
            matcher_t matcher(params);
            matcher.setCallback(matcher_t::callback_t::from<MatchingNode, &MatchingNode::publish>(this));

            // match the scans
            const double score = matcher.match(dst, src_, transform,
                                               use_odom ? initial_transform :
                                                          cslibs_math_2d::Transform2d::identity());

            std::cout << "matching took " << (cslibs_time::Time::now() - start).milliseconds() << "ms \n";
            std::cout << "score is " << score << "\n";
            std::cout << "transformation is " << transform << "\n";

            // publish results
            //if (score > 0.0) {

                // convert point clouds for publication
                sensor_msgs::PointCloud output, output_src, output_dst;
                output.header     = msg->header;
                output_src.header = msg->header;
                output_dst.header = msg->header;

                toPointCloud(*src_, output, transform.inverse());
                toPointCloud(*src_, output_src);
                toPointCloud(*dst,  output_dst);

                // publish transformed laserscan, src and dst
                pub_src_.publish(output_src);
                pub_dst_.publish(output_dst);
                pub_match_.publish(output);
            //}
        }

        // save the current laserscan as reference scan
        src_transform_ = dst_transform;
        src_           = dst;
    }

    void publish(
            const nav_msgs::OccupancyGrid::Ptr & msg)
    const
    {
        if (msg)
            pub_map_.publish(msg);
    }

    void toPointCloud(
            const cslibs_math_2d::Pointcloud2d & cloud,
            sensor_msgs::PointCloud            & output,
            const cslibs_math_2d::Transform2d  & transform = cslibs_math_2d::Transform2d::identity())
    const
    {
        for (const cslibs_math_2d::Point2d & point : cloud) {
            const cslibs_math_2d::Point2d p = transform * point;

            geometry_msgs::Point32 gp;
            gp.x = static_cast<float>(p(0));
            gp.y = static_cast<float>(p(1));
            output.points.emplace_back(gp);
        }
    }

    ros::NodeHandle                   nh_;
    ros::Subscriber                   sub_;
    ros::Publisher                    pub_src_;
    ros::Publisher                    pub_match_;
    ros::Publisher                    pub_dst_;
    ros::Publisher                    pub_map_;
    tf::TransformListener             tf_;

    std::string                       odom_frame_;
    cslibs_math_2d::Pointcloud2d::Ptr src_;
    tf::StampedTransform              src_transform_;
};
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "cslibs_ndt_2d_matching_node");
    cslibs_ndt_2d::MatchingNode node;
    ros::spin();

    return 0;
}
