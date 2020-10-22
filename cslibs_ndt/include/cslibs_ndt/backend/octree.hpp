#ifndef CSLIBS_NDT_BACKEND_OCTREE_HPP
#define CSLIBS_NDT_BACKEND_OCTREE_HPP

#include <cslibs_indexed_storage/backend/tags.hpp>
#include <cslibs_indexed_storage/backend/backend_traits.hpp>
#include <cslibs_indexed_storage/interface/data/data_interface.hpp>
#include <cslibs_indexed_storage/interface/data/align/aligned_allocator.hpp>

namespace cslibs_indexed_storage { namespace backend {
struct octree_tag {};
}}

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
namespace backend {

template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
class OcTree
{
public:
    using tag = cis::backend::octree_tag;

    using data_if = data_interface_t_;
    using data_storage_t = typename data_if::storage_type;
    using data_output_t = typename data_if::output_type;

    using index_if = index_interface_t_;
    using index_t = typename index_if::type;

    static constexpr auto on_duplicate_index_strategy =
            cis::option::get_option<cis::option::merge_strategy_opt, options_ts_...>::value;
    static constexpr std::size_t dimension = (1 << index_if::dimensions);

protected:
    struct Data
    {
        index_t index;
        data_storage_t data;
    };

    union Node {
        inline bool childExists(const unsigned int pos) const
        {
            assert (pos < dimension);
            return children_ && children_[pos];
        }

        inline bool hasChildren() const
        {
            if (!children_)
                return false;

            for (unsigned int i=0; i<dimension; ++i)
                if (children_[i])
                    return true;

            return false;
        }

        inline Node* createChild(const unsigned int pos)
        {
            assert (pos < dimension);
            // allocate children pointers
            if (!children_)
                allocateChildren();

            // create children
            assert (!children_[pos]);
            children_[pos] = new Node();

            return children_[pos];
        }

        inline Node* getChild(const unsigned int pos)
        {
            assert ((pos < dimension) && children_);
            assert (children_[pos]);
            return children_[pos];
        }

        inline const Node* getChild(const unsigned int pos) const
        {
            assert ((pos < dimension) && children_);
            assert (children_[pos]);
            return children_[pos];
        }

        inline void allocateChildren()
        {
            children_ = new Node*[dimension];
            for (unsigned int i=0; i<dimension; ++i)
                children_[i] = nullptr;
        }

        template<typename... Args>
        inline data_output_t& insert(const bool node_just_created, const index_t& index, Args&&... args)
        {
            // in new node, create new data, set index
            if (node_just_created) {
                data_ptr_ = new Data();

                auto& value = data_ptr_->data;
                value = data_if::create(std::forward<Args>(args)...);
                data_ptr_->index = index;
                return data_if::expose(value);
            }
            // in old node, merge data
            else {
                auto& value = data_ptr_->data;
                data_if::template merge<on_duplicate_index_strategy>(value, std::forward<Args>(args)...);
                return data_if::expose(value);
            }
        }

        inline data_storage_t* get()
        {
            return data_ptr_ ? &data_if::expose(data_ptr_->data) : nullptr;
        }

        template<typename Fn>
        inline void apply(const Fn& function)
        {
            if (data_ptr_)
                function(data_ptr_->index, data_if::expose(data_ptr_->data));
        }

        template<typename Fn>
        inline void apply(const Fn& function) const
        {
            if (data_ptr_)
                function(data_ptr_->index, data_if::expose(data_ptr_->data));
        }

        inline void deleteChildren()
        {
            if (children_) {
                for (unsigned int i=0; i<dimension; ++i) {
                    if (children_[i]) {
                        delete children_[i];
                    }
                }
                delete children_;
            }
        }

        inline void deleteData()
        {
            if (data_ptr_) {
                data_ptr_->index = index;
                data_if::deallocate(data_ptr_->data);
                delete data_ptr_;
            }
        }

        inline std::size_t byte_size() const
        {
            if (data_ptr_) {
                return sizeof(index_t) + data_if::byte_size(data_ptr_->data);
            }
            return 0;
        }

    private:
        Node** children_;
        Data*  data_ptr_ = nullptr;
    };

public:
    template<typename... Args>
    inline data_output_t& insert(const index_t& index, Args&&... args)
    {
        bool created_root = false;
        if (root_ == nullptr) {
            root_ = new Node();
            ++tree_size_;
            created_root = true;
        }

        return insert(root_, created_root, index, 0, args...);
    }

    inline data_output_t* get(const index_t& index)
    {
        if (root_)
            return get(root_, index, 0);

        return nullptr;
    }

    inline const data_output_t* get(const index_t& index) const
    {
        if (root_)
            return get(root_, index, 0);

        return nullptr;
    }

    template<typename Fn>
    inline void traverse(const Fn& function)
    {
        if (root_)
            traverse(root_, function, 0);
    }

    template<typename Fn>
    inline void traverse(const Fn& function) const
    {
        if (root_)
            traverse(root_, function, 0);
    }

    inline void clear()
    {
        if (root_) {
            clear(root_, 0);
            delete root_;
        }
        tree_size_ = 0;
    }

    virtual inline std::size_t byte_size() const
    {
        return root_ ? (sizeof(*this) + byte_size(root_,0,0)) : sizeof(*this);
    }

    inline std::size_t size() const
    {
        return tree_size_;
    }

private:
    template<typename... Args>
    inline data_output_t& insert(Node* node, const bool node_just_created,
                                 const index_t& index, const unsigned int depth, Args&&... args)
    {
        bool created_node = false;
        assert(node);

        // recurse down to last level
        if (depth < tree_depth_) {
            const unsigned int pos = computeChildIdx(index, tree_depth_ -1 - depth);
            if (!node->childExists(pos)) {
                node->createChild(pos);
                ++tree_size_;
                created_node = true;
            }
            return insert(node->getChild(pos), created_node, index, depth +1, args...);
        }

        // at last level, update node, end of recursion
        return node->insert(node_just_created, index, args...);
    }

    inline data_output_t* get(Node* node, const index_t& index, const unsigned int depth)
    {
        assert (node);

        // recurse down to last level
        if (depth < tree_depth_) {
            const unsigned int pos = computeChildIdx(index, tree_depth_ -1 - depth);
            if (!node->childExists(pos))
                return nullptr;

            return get(node->getChild(pos), index, depth +1);
        }

        // at last level, get data, end of recursion
        return node->get();
    }

    inline data_output_t* get(const Node* node, const index_t& index, const unsigned int depth) const
    {
        assert (node);

        // recurse down to last level
        if (depth < tree_depth_) {
            const unsigned int pos = computeChildIdx(index, tree_depth_ -1 - depth);
            if (!node->childExists(pos))
                return nullptr;

            return get(node->getChild(pos), index, depth +1);
        }

        // at last level, get data, end of recursion
        return node->get();
    }

    template<typename Fn>
    inline void traverse(Node* node, const Fn& function, const unsigned int depth)
    {
        assert (node);

        // recurse down to last level
        if (depth < tree_depth_) {
            for (std::size_t i = 0; i < dimension; ++i) {
                if (node->childExists(i))
                    traverse(node->getChild(i), function, depth + 1);
            }
        }

        // at last level, apply function, end of recursion
        else
            node->apply(function);
    }

    template<typename Fn>
    inline void traverse(const Node* node, const Fn& function, const unsigned int depth) const
    {
        assert (node);

        // recurse down to last level
        if (depth < tree_depth_) {
            for (std::size_t i = 0; i < dimension; ++i) {
                if (node->childExists(i))
                    traverse(node->getChild(i), function, depth + 1);
            }
        }

        // at last level, apply function, end of recursion
        else
            node->apply(function);
    }

    template<typename Fn>
    inline void clear(Node* node, const unsigned int depth)
    {
        assert (node);

        // recurse down to last level
        if (depth < tree_depth_) {
            for (std::size_t i = 0; i < dimension; ++i) {
                if (node->childExists(i)) {
                    clear(node->getChild(i), depth + 1);
                    node->deleteChildren();
                }
            }
        }

        // at last level, clear data, end of recursion
        else {
            node->clear();
            node->deleteData();
        }
    }

    inline std::size_t byte_size(const Node* node, const unsigned int depth, std::size_t bytes) const
    {
        assert (node);
        bytes += sizeof(Node);

        // recurse down to last level
        if (depth < tree_depth_) {
            if (node->hasChildren()) {
                bytes += dimension * sizeof(Node*);

                for (std::size_t i = 0; i < dimension; ++i) {
                    if (node->childExists(i)) {
                        bytes += byte_size(node->getChild(i), depth +1, 0);
                    }
                }
            }
            return bytes;
        }

        // at last level, get byte size of data, end of recursion
        return bytes + node->byte_size();
    }

    inline unsigned int computeChildIdx(const index_t& index, const unsigned int depth)
    {
        unsigned int pos = 0;
        for (std::size_t i = 0; i < index_if::dimensions; ++i) {
            const int key = index_if::access(i, index) + tree_max_val_;
            if (key & (1 << depth))
                pos += (1 << i);
        }

        return pos;
    }

protected:
    Node*              root_ = nullptr;
    std::size_t        tree_size_;
    const unsigned int tree_depth_ = 16;         //TODO???
    const int          tree_max_val_ = 32768;    // = 2^15 = (1 << (tree_depth_ - 1))
};

}
}


#endif // CSLIBS_NDT_BACKEND_OCTREE_HPP
