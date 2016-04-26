#ifndef NDT_KDTREE_HPP
#define NDT_KDTREE_HPP

#include <cstddef>

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
    typedef typename std::vector<DistributionType>  DistributionSetType;
    typedef typename DistributionType::MatrixType   CovarianceMatrixType;

    typedef typename DistributionType::PointType    PointType;

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
        if(!distribution || distribution.getN() != points.size()) {
            distribution.reset(new DistributionType);
            for(const PointType &p : points) {
                distribution.add(p);
            }
        }
        return distribution.get();
    }
};

template<std::size_t Dim>
struct KDTreeIndex
{
    typedef std::array<std::size_t, Dim> IndexType;
    typedef std::array<double, Dim>      ResolutionType;
    typedef KDTreeNode<Dim>              NodeType;

    KDTreeIndex()
    {
        resolution.fill(0.5);
    }

    KDTreeIndex(const ResolutionType &_resolution) :
        resolution(_resolution)
    {
    }

    typename NodeType::Ptr create(const typename NodeType::PointType &_p)
    {
        NodeType *node = new NodeType;
        /// ... extra material
        return NodeType::Ptr(node);
    }

    ResolutionType resolution;

};


}
}
#endif // NDT_KDTREE_HPP
