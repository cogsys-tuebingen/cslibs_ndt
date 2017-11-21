#ifndef VIS_KDTREE_HPP
#define VIS_KDTREE_HPP

/// SYSTEM
#include <opencv2/opencv.hpp>

/// PROJECT
#include <ndt/tree/kdtree.hpp>

namespace ndt {
namespace visualization {
typedef Eigen::Vector2d                                     Point2D;
typedef ndt::math::Distribution<2, true>                    Distribution2D;
typedef ndt::tree::Index<2>                                 KDIndex2D;
typedef ndt::tree::NodeData<2>                              KDNodeData2D;
typedef kdtree::unbuffered::KDTree<KDIndex2D, KDNodeData2D> KDTree2D;
typedef KDTree2D::NodeType                                  KDNode2D;


void renderTree(KDTree2D::Ptr   &tree,
                KDIndex2D         &index,
                const Point2D     &min,
                const Point2D     &max,
                cv::Mat           &dst)
{
    double scale_x = fabs((max - min)(0)) / dst.cols;
    double scale_y = fabs((max - min)(1)) / dst.rows;

    cv::Mat samples(dst.rows, dst.cols, CV_64FC1, cv::Scalar());
    double max_value = std::numeric_limits<double>::lowest();
//#pragma omp parallel for
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            Point2D p = min + Point2D(scale_x * j, scale_y * i);


            KDNode2D *node = tree->find(index.create(p));
            if(node == nullptr)
                continue;
            if(!node->data.distribution)
                continue;

            Distribution2D &distribution = *(node->data.distribution);

            if(distribution.getN() < 3)
                continue;

            double value = distribution.sampleNonNormalized(p);
            if(value > max_value) {
                max_value = value;
            }
            samples.at<double>(dst.rows - 1 - i,j) = value;
        }
    }

#pragma omp parallel for
    for(int i = 0 ; i < dst.rows ; ++i) {
        for(int j = 0 ; j < dst.cols ; ++j) {
            cv::Vec3b &pix = dst.at<cv::Vec3b>(i,j);
            pix[0] = samples.at<double>(i,j) / max_value * 255;
            pix[1] = samples.at<double>(i,j) / max_value * 255;
            pix[2] = samples.at<double>(i,j) / max_value * 255;
        }
    }
}

//void renderClusterDistributions(KDTreeType::Ptr   &tree,
//                                IndexType         &index,
//                                const Point2D     &min,
//                                const Point2D     &max,
//                                cv::Mat           &dst)
//{
//    double scale_x = fabs((max - min)(0)) / dst.cols;
//    double scale_y = fabs((max - min)(1)) / dst.rows;

//    KDTreeInterface2D::DistributionMapType distributions;
//    interface.getClusterDistributions(tree, distributions);

//    cv::Mat samples(dst.rows, dst.cols, CV_64FC1, cv::Scalar());
//    double max_value = std::numeric_limits<double>::lowest();
//    for(int i = 0 ; i < dst.rows ; ++i) {
//        for(int j = 0 ; j < dst.cols ; ++j) {
//            Point2D p = min + Point2D(scale_x * j, scale_y * i);

//            double value = 0.0;
//            for(auto distributions_entry : distributions) {
//                KDTreeInterface2D::DistributionType &distribution = distributions_entry.second;
//                if(distribution.getN() < 3)
//                    continue;
//                    value += distribution.sampleNonNormalized(p);
//            }

//            if(value > max_value) {
//                max_value = value;
//            }
//            samples.at<double>(dst.rows - 1 - i,j) = value;
//        }
//    }

//#pragma omp parallel for
//    for(int i = 0 ; i < dst.rows ; ++i) {
//        for(int j = 0 ; j < dst.cols ; ++j) {
//            cv::Vec3b &pix = dst.at<cv::Vec3b>(i,j);
//            pix[0] = samples.at<double>(i,j) / max_value * 255;
//            pix[1] = samples.at<double>(i,j) / max_value * 255;
//            pix[2] = samples.at<double>(i,j) / max_value * 255;
//        }
//    }
//}

}
}

#endif // VIS_KDTREE_HPP
