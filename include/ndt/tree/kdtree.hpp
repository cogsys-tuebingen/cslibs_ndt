#ifndef NDT_KDTREE_HPP
#define NDT_KDTREE_HPP

#include <cstddef>

#include <ndt/data/pointcloud.hpp>
#include <ndt/math/distribution.hpp>
#include <kdtree/kdtree.hpp>


namespace ndt {
namespace tree {
template<std::size_t Dim>
struct KDTreeNode : public kdtree::KDTreeNode<int, Dim>
{
    typedef kdtree::KDTree<int, Dim>                KDTreeType;
    typedef kdtree::KDTreeNode<int, Dim>            NodeBase;
    typedef math::Distribution<Dim, true>           DistributionType;
    typedef typename DistributionType::MatrixType   CovarianceMatrixType;
    typedef data::Pointcloud<Dim>                   PointCloudType;
    typedef typename PointCloudType::PointType      PointType;

    std::vector<PointType>          points;
    typename DistributionType::Ptr  distribution;

    KDTreeNode() :
        NodeBase(),
        distribution(nullptr)
    {
    }

    KDTreeNode(const KDTreeNode &other) :
        NodeBase(other),
        points(other.points),
        distribution(other.distribution)
    {
    }

    inline typename NodeBase::Ptr clone() const override
    {
        typename NodeBase::Ptr node(new KDTreeNode<Dim>(*this));
        node->left  = NodeBase::left;
        node->right = NodeBase::right;
        return node;
    }

    inline typename NodeBase::Ptr copy() const override
    {
        return typename NodeBase::Ptr(new KDTreeNode(*this));
    }

    inline void overwrite(const typename NodeBase::Ptr &other) override
    {
        KDTreeNode *other_ptr = (KDTreeNode*) other.get();
        points.insert(points.end(),
                      other_ptr->points.begin(),
                      other_ptr->points.end());
    }

    inline void split(const typename NodeBase::Ptr &other) override
    {
        NodeBase::split(other);
        points.clear();
    }

    inline DistributionType* get()
    {
        if(!distribution || distribution->getN() != points.size()) {
            distribution.reset(new DistributionType);
            for(const PointType &p : points) {
                distribution->add(p);
            }
        }
        return distribution.get();
    }
};

template<std::size_t Dim>
struct KDTreeInterface
{
    typedef std::array<int,    Dim>                 IndexType;
    typedef std::array<double, Dim>                 ResolutionType;
    typedef KDTreeNode<Dim>                         NodeType;
    typedef typename NodeType::KDTreeType           KDtreeType;
    typedef data::Pointcloud<Dim>                   PointCloudType;
    typedef typename PointCloudType::PointType      PointType;
    typedef math::Distribution<Dim, true>           DistributionType;

    KDTreeInterface()
    {
        resolution.fill(0.5);
    }

    KDTreeInterface(const ResolutionType &_resolution) :
        resolution(_resolution)
    {
    }

    inline typename NodeType::Ptr create(const typename NodeType::PointType &_p)
    {
        NodeType *node = new NodeType;
        node->points.emplace_back(_p);
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            node->index[i] = floor(_p(i) / resolution[i]);
        }
        return typename NodeType::Ptr(node);
    }

    inline DistributionType * get(const typename NodeType::PointType &_p,
                                  typename KDtreeType::Ptr &_tree)
    {
        if(!_tree)
            return nullptr;

        IndexType index;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            index[i] = floor(_p(i) / resolution[i]);
        }
        typename NodeType::Ptr node;
        if(_tree->find(index, node)) {
            return ((NodeType *) node.get())->get();
        }

        return nullptr;
    }

    inline void insert(const typename PointCloudType::Ptr &_point_cloud,
                       typename KDtreeType::Ptr &_tree)
    {
        const PointCloudType &point_cloud = *_point_cloud;
        insert(point_cloud, _tree);
    }

    inline void insert(const PointCloudType &_point_cloud,
                       typename KDtreeType::Ptr &_tree)
    {
        _tree.reset(new KDtreeType);
        for(std::size_t i = 0 ; i < _point_cloud.size ; ++i) {
            if(_point_cloud.mask[i] == PointCloudType::VALID)
                _tree->insertNode(create(_point_cloud.points[i]));
        }
    }

    ResolutionType resolution;

};


}
}
#endif // NDT_KDTREE_HPP
