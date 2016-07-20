#ifndef NDT_KDTREE_HPP
#define NDT_KDTREE_HPP

#include <cstddef>

#include <ndt/data/pointcloud.hpp>
#include <ndt/math/distribution.hpp>
#include <kdtree/kdtree.hpp>
#include <kdtree/kdtree_clustering.hpp>
#include <mutex>

namespace ndt {
namespace tree {
template<std::size_t Dim>
struct Index {
    using PointType = Eigen::Matrix<double, Dim, 1>;
    using BinType   = std::array<double, Dim>;
    using Type = std::array<int,Dim>;
    using PivotType = double;

    static constexpr std::size_t Dimension = Dim;

    const BinType sizes;

    Index() = default;

    Index(const BinType &sizes) :
        sizes(sizes)
    {
    }

    static Type fromPoint(const PointType &_p,
                          const BinType &_sizes)
    {
        Type res;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            res[i] = static_cast<int>(floor(_p(i) / _sizes[i]));
        }
        return res;
    }

    inline constexpr Type create(const PointType &_p) const
    {
        return fromPoint(_p, sizes);
    }
};

template<std::size_t Dim>
struct NodeData {
    typedef math::Distribution<Dim, true> DistributionType;
    typename DistributionType::Ptr   distribution;

    NodeData() :
        distribution(new DistributionType)
    {
    }

    inline NodeData(const typename Index<Dim>::PointType &_p) :
        distribution(new DistributionType)
    {
        distribution->add(_p);
    }

    inline void merge(const NodeData &other)
    {
        if(other.distribution->getN() > 0)
            *distribution += (*other.distribution);
    }
};

//template<std::size_t Dim>
//struct KDTreeInterface
//{
//    typedef std::array<int,    Dim>                   IndexType;
//    typedef std::array<double, Dim>                   ResolutionType;
//    typedef KDTreeNode<Dim>                           NodeType;
//    typedef typename NodeType::KDTreeType             KDtreeType;
//    typedef typename NodeType::KDTreeClusteringType   KDTreeClusteringType;
//    typedef data::Pointcloud<Dim>                     PointCloudType;
//    typedef typename PointCloudType::PointType        PointType;
//    typedef math::Distribution<Dim, true>             DistributionType;
//    typedef typename std::map<int, DistributionType>  DistributionMapType;

//    KDTreeInterface()
//    {
//        resolution.fill(0.5);
//    }

//    KDTreeInterface(const ResolutionType &_resolution) :
//        resolution(_resolution)
//    {
//    }

//    inline typename NodeType::Ptr create(const typename NodeType::PointType &_p)
//    {
//        NodeType *node = new NodeType(_p);

//        for(std::size_t i = 0 ; i < Dim ; ++i) {
//            node->index[i] = floor(_p(i) / resolution[i]);
//        }

//        return typename NodeType::Ptr(node);
//    }

//    inline DistributionType * get(const typename NodeType::PointType &_p,
//                                  typename KDtreeType::Ptr &_tree)
//    {
//        if(!_tree)
//            return nullptr;

//        IndexType index;
//        for(std::size_t i = 0 ; i < Dim ; ++i) {
//            index[i] = floor(_p(i) / resolution[i]);
//        }

//        typename NodeType::Ptr node;
//        if(_tree->find(index, node)) {
//            return ((NodeType*) node.get())->distribution.get();
//        }

//        return nullptr;
//    }

//    inline void insert(const typename PointCloudType::Ptr &_point_cloud,
//                       typename KDtreeType::Ptr &_tree)
//    {
//        const PointCloudType &point_cloud = *_point_cloud;
//        insert(point_cloud, _tree);
//    }

//    inline void insert(const PointCloudType &_point_cloud,
//                       typename KDtreeType::Ptr &_tree)
//    {
//        _tree.reset(new KDtreeType);
//        for(std::size_t i = 0 ; i < _point_cloud.size ; ++i) {
//            if(_point_cloud.mask[i] == PointCloudType::VALID)
//                _tree->insertNode(create(_point_cloud.points[i]));
//        }
//    }

//    inline void cluster(typename KDtreeType::Ptr &_tree)
//    {
//        if(!_tree)
//            return;
//        KDTreeClusteringType clustering(_tree);
//        clustering.cluster();
//    }

//    inline void getClusterDistributions(const typename KDtreeType::Ptr &_tree,
//                                        DistributionMapType &_distributions)
//    {
//        std::vector<typename NodeType::Ptr> leaves;
//        _tree->getLeaves(leaves);
//        for(typename NodeType::Ptr &leaf : leaves) {
//            NodeType *node = leaf.get();
//            DistributionType &distr = _distributions[node->cluster];
//            if(node->distribution)
//                distr += *(node->distribution);
//        }
//    }


//    ResolutionType resolution;

//};


}
}
#endif // NDT_KDTREE_HPP
