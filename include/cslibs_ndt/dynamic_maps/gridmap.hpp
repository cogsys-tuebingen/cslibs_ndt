#ifndef CSLIBS_GRIDMAPS_DYNAMIC_GRIDMAP_HPP
#define CSLIBS_GRIDMAPS_DYNAMIC_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_2d/linear/pose.hpp>
#include <cslibs_math_2d/linear/point.hpp>

#include <cslibs_gridmaps/dynamic_maps/algorithms/bresenham.hpp>
#include <cslibs_gridmaps/dynamic_maps/chunk.hpp>

#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/mod.hpp>
#include <cslibs_math/common/div.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_gridmaps {
namespace dynamic_maps {
template<typename T>
class Gridmap
{
public:
    using Ptr                   = std::shared_ptr<Gridmap<T>>;
    using gridmap_t             = Gridmap<T>;
    using pose_t                = cslibs_math_2d::Pose2d;
    using index_t               = std::array<int, 2>;
    using mutex_t               = std::mutex;
    using lock_t                = std::unique_lock<mutex_t>;
    using chunk_t               = Chunk<T>;
    using storage_t             = cis::Storage<chunk_t, index_t, cis::backend::kdtree::KDTree>;
    using line_iterator_t       = algorithms::Bresenham<T>;
    using const_line_iterator_t = algorithms::Bresenham<T const>;
    using get_chunk_t           = delegate<chunk_t*(const index_t &)>;

    Gridmap(const pose_t        &origin,
            const double         resolution,
            const double         chunk_resolution,
            const T             &default_value) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        chunk_size_(static_cast<int>(chunk_resolution * resolution_inv_)),
        default_value_(default_value),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_chunk_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_chunk_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        storage_(new storage_t),
        height_(chunk_size_),
        width_(chunk_size_)
    {
    }

    Gridmap(const double origin_x,
            const double origin_y,
            const double origin_phi,
            const double resolution,
            const double chunk_resolution,
            const T &default_value) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        chunk_size_(static_cast<int>(chunk_resolution * resolution_inv_)),
        default_value_(default_value),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        min_chunk_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_chunk_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        storage_(new storage_t),
        height_(chunk_size_),
        width_(chunk_size_)
    {
    }

    inline cslibs_math_2d::Point2d getMin() const
    {
        return fromIndex({0,0});
    }

    inline cslibs_math_2d::Point2d getMax() const
    {
        return fromIndex({static_cast<int>(width_ - 1),
                          static_cast<int>(height_ - 1)});
    }

    inline cslibs_math_2d::Pose2d getOrigin() const
    {
        cslibs_math_2d::Transform2d origin = w_T_m_;
        origin.translation() = fromIndex({{0,0}});
        return origin;
    }

    inline cslibs_math_2d::Pose2d getInitialOrigin() const
    {
        return w_T_m_;
    }

    /// build a handle
    inline T& at(const std::size_t idx, const std::size_t idy)
    {
        if(idx >= width_ || idy >= height_) {
            throw std::runtime_error("[GridMap] : Invalid Index!");
        }

        const index_t index             = {idx + min_index_[0], idy + min_index_[1]};
        const index_t chunk_index       = toChunkIndex(index);
        const index_t local_chunk_index = toLocalChunkIndex(index);

        return getAllocateChunk(chunk_index)->at(local_chunk_index);
    }

    inline T at(const std::size_t idx, const std::size_t idy) const
    {
        if(idx >= width_ || idy >= height_) {
            throw std::runtime_error("[GridMap] : Invalid Index!");
        }

        const index_t index             = {static_cast<int>(idx) + min_index_[0],
                                           static_cast<int>(idy) + min_index_[1]};
        const index_t chunk_index       = toChunkIndex(index);
        const index_t local_chunk_index = toLocalChunkIndex(index);

        const chunk_t *chunk = lookupChunk(chunk_index);
        return chunk != nullptr ? chunk->at(local_chunk_index) : default_value_;
    }

    inline T& at(const cslibs_math_2d::Point2d &point)
    {
        const index_t index             = toIndex(point);
        const index_t chunk_index       = toChunkIndex(index);
        const index_t local_chunk_index = toLocalChunkIndex(index);

        return getAllocateChunk(chunk_index)->at(local_chunk_index);
    }

    inline T at(const cslibs_math_2d::Point2d &point) const
    {
        const index_t index             = toIndex(point);
        const index_t chunk_index       = toChunkIndex(index);
        const index_t local_chunk_index = toLocalChunkIndex(index);
        const chunk_t *chunk = lookupChunk(chunk_index);
        return chunk != nullptr ? chunk->at(local_chunk_index) : default_value_;
    }

    inline line_iterator_t getLineIterator(const index_t &start_index,
                                           const index_t &end_index)
    {
        return line_iterator_t(start_index,
                               end_index,
                               chunk_size_,
                               default_value_,
                               get_chunk_t::template from<gridmap_t, &gridmap_t::getAllocateChunk>(this));
    }

    inline line_iterator_t getLineIterator(const cslibs_math_2d::Point2d &start,
                                           const cslibs_math_2d::Point2d &end)
    {

        const index_t start_index = toIndex(start);
        const index_t end_index   = toIndex(end);

        return  line_iterator_t(start_index, end_index,
                                chunk_size_,
                                default_value_,
                                get_chunk_t::template from<gridmap_t, &gridmap_t::getAllocateChunk>(this));
    }

    inline line_iterator_t getLineIterator(const index_t &start_index,
                                           const index_t &end_index) const
    {
        return line_iterator_t(start_index,
                               end_index,
                               chunk_size_,
                               default_value_,
                               get_chunk_t::template from<gridmap_t, &gridmap_t::getAllocateChunk>(this));
    }

    inline line_iterator_t getLineIterator(const cslibs_math_2d::Point2d &start,
                                           const cslibs_math_2d::Point2d &end) const
    {

        const index_t start_index = toIndex(start);
        const index_t end_index   = toIndex(end);
        return  line_iterator_t(start_index, end_index,
                                chunk_size_,
                                default_value_,
                                get_chunk_t::template from<gridmap_t, &gridmap_t::getAllocateChunk>(this));
    }


    inline index_t getMinChunkIndex() const
    {
        return min_chunk_index_;
    }

    inline index_t getMaxChunkIndex() const
    {
        return max_chunk_index_;
    }

    inline chunk_t const * lookupChunk(const index_t &chunk_index) const
    {
        lock_t l(storage_mutex_);
        return storage_->get(chunk_index);
    }

    inline chunk_t* getChunk(const index_t &chunk_index) const
    {
        lock_t l(storage_mutex_);
        chunk_t *chunk = storage_->get(chunk_index);
        return chunk;
    }

    inline chunk_t* getAllocateChunk(const index_t &chunk_index) const
    {
        lock_t l(storage_mutex_);
        chunk_t *chunk = storage_->get(chunk_index);
        if(chunk == nullptr) {
            chunk = &(storage_->insert(chunk_index, chunk_t(chunk_size_, default_value_)));
        }
        updateChunkIndices(chunk_index);
        return chunk;
    }

    inline double getResolution() const
    {
        return resolution_;
    }

    inline std::size_t getChunkSize() const
    {
        return static_cast<std::size_t>(chunk_size_);
    }

    inline std::size_t getHeight() const
    {
        return height_;
    }

    inline std::size_t getWidth() const
    {
        return width_;
    }

    inline index_t getMaxIndex() const
    {
        return {(max_chunk_index_[0] - min_chunk_index_[0] + 1) * chunk_size_ - 1,
                (max_chunk_index_[1] - min_chunk_index_[1] + 1) * chunk_size_ - 1};
    }


protected:
    const double                      resolution_;
    const double                      resolution_inv_;
    const int                         chunk_size_;
    const T                           default_value_;
    cslibs_math_2d::Transform2d       w_T_m_;
    cslibs_math_2d::Transform2d       m_T_w_;


    mutable index_t                    min_chunk_index_;
    mutable index_t                    max_chunk_index_;
    mutable index_t                    min_index_;
    mutable mutex_t                    storage_mutex_;
    mutable std::shared_ptr<storage_t> storage_;
    mutable std::size_t                height_;
    mutable std::size_t                width_;

    inline void updateChunkIndices(const index_t &chunk_index) const
    {
        min_chunk_index_    = cslibs_math::common::min(min_chunk_index_, chunk_index);
        max_chunk_index_    = cslibs_math::common::max(max_chunk_index_, chunk_index);
        min_index_          = min_chunk_index_ * chunk_size_;
        width_  = (max_chunk_index_[0] - min_chunk_index_[0] + 1) * chunk_size_;
        height_ = (max_chunk_index_[1] - min_chunk_index_[1] + 1) * chunk_size_;
    }

    inline index_t toChunkIndex(const index_t &index) const
    {
        return {{cslibs_math::common::div(index[0], chunk_size_),
                 cslibs_math::common::div(index[1], chunk_size_)}};
    }

    inline index_t toLocalChunkIndex(const index_t &index) const
    {
        return {{cslibs_math::common::mod(index[0], chunk_size_),
                 cslibs_math::common::mod(index[1], chunk_size_)}};
    }

    inline index_t toIndex(const cslibs_math_2d::Point2d &p_w) const
    {
        /// offset and rounding correction!
        const cslibs_math_2d::Point2d p_m = m_T_w_ * p_w;
        return {{static_cast<int>(p_m(0) * resolution_inv_),
                 static_cast<int>(p_m(1) * resolution_inv_)}};
    }

    inline cslibs_math_2d::Point2d fromIndex(const index_t &i) const
    {
        return w_T_m_ * cslibs_math_2d::Point2d((i[0] + min_index_[0])  * resolution_,
                                                (i[1] + min_index_[0])  * resolution_);
    }

};
}
}



#endif // CSLIBS_GRIDMAPS_DYNAMIC_GRIDMAP_HPP
