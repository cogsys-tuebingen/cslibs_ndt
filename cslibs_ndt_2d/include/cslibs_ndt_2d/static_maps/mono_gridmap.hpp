#ifndef CSLIBS_NDT_2D_STATIC_MAPS_FLAT_GRIDMAP_HPP
#define CSLIBS_NDT_2D_STATIC_MAPS_FLAT_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_2d/linear/pose.hpp>
#include <cslibs_math_2d/linear/point.hpp>

#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/linear/pointcloud.hpp>
#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_2d {
namespace static_maps {
namespace mono {
template <typename T>
class EIGEN_ALIGN16 Gridmap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Gridmap<T>>;

    using ConstPtr                          = std::shared_ptr<const Gridmap<T>>;
    using Ptr                               = std::shared_ptr<Gridmap<T>>;
    using pose_t                            = cslibs_math_2d::Pose2d<T>;
    using transform_t                       = cslibs_math_2d::Transform2d<T>;
    using point_t                           = cslibs_math_2d::Point2d<T>;
    using index_t                           = std::array<int, 2>;
    using size_t                            = std::array<std::size_t, 2>;
    using size_m_t                          = std::array<T, 2>;
    using distribution_t                    = cslibs_ndt::Distribution<T,2>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_const_ptr_t  = std::shared_ptr<distribution_storage_t const>;

    inline Gridmap(const pose_t  &origin,
                   const T       &resolution,
                   const size_t  &size,
                   const index_t &min_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
        (size[1] + 1) * resolution}},
        min_index_(min_index),
        max_index_{{min_index[0] + static_cast<int>(size[0]),
        min_index[1] + static_cast<int>(size[1])}},
        storage_(distribution_storage_ptr_t(new distribution_storage_t))
    {
        storage_->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_->template set<cis::option::tags::array_offset>(min_index[0], min_index[1]);
    }

    inline Gridmap(const T &origin_x,
                   const T &origin_y,
                   const T &origin_phi,
                   const T &resolution,
                   const size_t &size,
                   const index_t &min_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
        (size[1] + 1) * resolution}},
        min_index_(min_index),
        max_index_{{min_index[0] + static_cast<int>(size[0]),
        min_index[1] + static_cast<int>(size[1])}},
        storage_(distribution_storage_ptr_t(new distribution_storage_t))
    {
        storage_->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_->template set<cis::option::tags::array_offset>(min_index[0], min_index[1]);
    }

    inline Gridmap(const Gridmap &other) :
        resolution_(other.resolution_),
        resolution_inv_(other.resolution_inv_),
        w_T_m_(other.w_T_m_),
        m_T_w_(other.m_T_w_),
        size_(other.size_),
        size_m_(other.size_m_),
        min_index_(other.min_index_),
        max_index_(other.max_index_),
        storage_(distribution_storage_ptr_t(new distribution_storage_t(*other.storage_)))
    {
    }

    inline Gridmap(Gridmap &&other) :
        resolution_(other.resolution_),
        resolution_inv_(other.resolution_inv_),
        w_T_m_(std::move(other.w_T_m_)),
        m_T_w_(std::move(other.m_T_w_)),
        size_(other.size_),
        size_m_(other.size_m_),
        min_index_(other.min_index_),
        max_index_(other.max_index_),
        storage_(other.storage_)
    {
    }

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    inline point_t getMin() const
    {
        return point_t(min_index_[0] * resolution_,
                min_index_[1] * resolution_);
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    inline point_t getMax() const
    {
        return point_t((max_index_[0] + 1) * resolution_,
                (max_index_[1] + 1) * resolution_);
    }

    /**
     * @brief Get the origin.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        transform_t origin = w_T_m_;
        origin.translation() += getMin();
        return origin;
    }

    /**
     * @brief Get the initial origin of the map.
     * @return the inital origin
     */
    inline pose_t getInitialOrigin() const
    {
        return w_T_m_;
    }

    inline index_t getMinIndex() const
    {
        return min_index_;
    }

    inline index_t getMaxIndex() const
    {
        return max_index_;
    }

    inline void insert(const point_t &p)
    {
        index_t i;
        if(!toIndex(p, i))
            return;

        distribution_t *distribution = getAllocate(i);
        distribution->data().add(p);
    }

    inline T sample(const point_t &p) const
    {
        index_t i;
        return toIndex(p, i) ? sample(p, i) : 0.0;
    }

    inline T sample(const point_t &p,
                    const index_t &i) const
    {
        distribution_t *distribution;
        distribution = storage_->get(i);

        return distribution ? distribution->data().sample(p) : 0.0;
    }

    inline T sampleNonNormalized(const point_t &p) const
    {
        index_t i;
        return toIndex(p, i) ? sampleNonNormalized(p, i) : 0.0;
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const index_t &i) const
    {
        distribution_t *distribution  = storage_->get(i);
        return distribution ? distribution->data().sampleNonNormalized(p) : 0.0;
    }

    inline distribution_t* get(const point_t &p) const
    {
        index_t i;
        if(!toIndex(p,i))
            return nullptr;

        return storage_->get(i);
    }

    inline const distribution_t* getDistribution(const index_t &i) const
    {
        return getAllocate(i);
    }

    inline distribution_t* getDistribution(const index_t &i)
    {
        return getAllocate(i);
    }

    inline T getResolution() const
    {
        return resolution_;
    }

    inline T getHeight() const
    {
        return size_[1] * resolution_;
    }

    inline T getWidth() const
    {
        return size_[0] * resolution_;
    }

    inline size_t getSize() const
    {
        return size_;
    }

    inline index_t getMinBundleIndex() const
    {
        return min_index_;
    }

    template <typename Fn>
    inline void traverse(const Fn& function) const
    {
        return storage_->traverse(function);
    }

    inline void getIndices(std::vector<index_t> &indices) const
    {
        auto add_index = [&indices](const index_t &i, const distribution_t &) {
            indices.emplace_back(i);
        };
        storage_->traverse(add_index);
    }

    inline std::size_t getByteSize() const
    {
        return sizeof(*this) +
                storage_->byte_size();
    }

    inline virtual bool validate(const pose_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w.translation();
        const index_t index = {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                                static_cast<int>(std::floor(p_m(1) * resolution_inv_))}};
        return (index[0] >= min_index_[0] && index[0] <= max_index_[0] ) &&
               (index[1] >= min_index_[1] && index[1] <= max_index_[1] );
    }

protected:
    const T                                    resolution_;
    const T                                    resolution_inv_;
    const transform_t                          w_T_m_;
    const transform_t                          m_T_w_;
    const size_t                               size_;
    const size_m_t                             size_m_;
    const index_t                              min_index_;
    const index_t                              max_index_;

    mutable distribution_storage_ptr_t         storage_;

    inline distribution_t *getAllocate(const index_t &i) const
    {
        distribution_t *distribution = storage_->get(i);

        auto allocate = [this, &i]() {
            return &(storage_->insert(i, distribution_t()));
        };
        return distribution ? distribution : allocate();
    }

    inline bool toIndex(const point_t &p_w,
                        index_t &index) const
    {
        const point_t p_m = m_T_w_ * p_w;
        index = {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                  static_cast<int>(std::floor(p_m(1) * resolution_inv_))}};
        return (index[0] >= min_index_[0] && index[0] <= max_index_[0] ) &&
               (index[1] >= min_index_[1] && index[1] <= max_index_[1] );
    }

    inline void fromIndex(const index_t &i,
                          point_t &p_w) const
    {
        p_w = w_T_m_ * point_t(i[0] * resolution_,
                               i[1] * resolution_);
    }
};
}
}
}

#endif // CSLIBS_NDT_2D_STATIC_MAPS_FLAT_GRIDMAP_HPP
