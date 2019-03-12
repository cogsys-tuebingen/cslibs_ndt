#ifndef CSLIBS_NDT_MAP_GENERIC_MAP_HPP
#define CSLIBS_NDT_MAP_GENERIC_MAP_HPP

#include <cslibs_ndt/map/abstract_map.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = tags::default_types<option_t>::template default_backend_t>
class EIGEN_ALIGN16 GenericMap : public AbstractMap<option_t, Dim, data_t, T, backend_t> {};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 GenericMap<tags::static_map, Dim, data_t, T, backend_t> : public AbstractMap<tags::static_map, Dim, data_t, T, backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t   = Eigen::aligned_allocator<GenericMap<tags::static_map, Dim, data_t, T, backend_t>>;

    using ConstPtr      = std::shared_ptr<const GenericMap<tags::static_map, Dim, data_t, T, backend_t>>;
    using Ptr           = std::shared_ptr<GenericMap<tags::static_map, Dim, data_t, T, backend_t>>;

    using base_t = AbstractMap<tags::static_map, Dim, data_t, T, backend_t>;
    using typename base_t::pose_t;
    using typename base_t::transform_t;
    using typename base_t::point_t;
    using typename base_t::pointcloud_t;
    using typename base_t::index_t;
    using typename base_t::index_list_t;
    using typename base_t::distribution_t;
    using typename base_t::distribution_storage_t;
    using typename base_t::distribution_storage_ptr_t;
    using typename base_t::distribution_storage_array_t;
    using typename base_t::distribution_bundle_t;
    using typename base_t::distribution_const_bundle_t;
    using typename base_t::distribution_bundle_storage_t;
    using typename base_t::distribution_bundle_storage_ptr_t;

    using size_t        = std::array<std::size_t,Dim>;
    using size_m_t      = std::array<T,Dim>;

    inline GenericMap(const pose_t  &origin,
                      const T       &resolution,
                      const size_t  &size,
                      const index_t &min_bundle_index) :
        base_t(origin, resolution,
               min_bundle_index,
               (min_bundle_index + size * 2 - 1)),
        size_(size),
        size_m_((size + 1) * resolution)
    {
        index_t offset;
        for (std::size_t i=0; i<Dim; ++i)
            offset[i] = cslibs_math::common::div<int>(min_bundle_index[i], 2);
        this->storage_[0]->template set<cis::option::tags::array_size>(std::forward(size));
        this->storage_[0]->template set<cis::option::tags::array_offset>(std::forward(offset));

        for(std::size_t i=1 ; i<Dim; ++i) {
            this->storage_[i]->template set<cis::option::tags::array_size>(std::forward(size+1));
            this->storage_[i]->template set<cis::option::tags::array_offset>(offset);
        }

        this->bundle_storage_->template set<cis::option::tags::array_size>(std::forward(size*2));
        this->bundle_storage_->template set<cis::option::tags::array_offset>(std::forward(min_bundle_index));
    }

    inline GenericMap(const pose_t &origin,
                      const T      &resolution,
                      const size_t &size,
                      const distribution_bundle_storage_ptr_t &bundles,
                      const distribution_storage_array_t      &storage,
                      const index_t                           &min_bundle_index) :
        base_t(origin, resolution,
               min_bundle_index,
               (min_bundle_index + size * 2 - 1),
               bundles, storage),
        size_(size),
        size_m_((size + 1) * resolution)
    {
    }

    inline GenericMap(const GenericMap &other) :
        base_t(other),
        size_(other.size_),
        size_m_(other.size_m_)
    {
    }

    inline GenericMap(GenericMap &&other) :
        base_t(other),
        size_(other.size_),
        size_m_(other.size_m_)
    {
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return valid(bi) ? this->getAllocate(bi) : nullptr;
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return valid(bi) ? this->getAllocate(bi) : nullptr;
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        index_t bi;
        if (!this->toBundleIndex(p, bi))
            return nullptr;

        return this->getAllocate(bi);
    }

    inline const distribution_bundle_t* get(const point_t &p) const
    {
        index_t bi;
        if (!this->toBundleIndex(p, bi))
            return nullptr;

        return this->bundle_storage_->get(bi);
    }

    inline const distribution_bundle_t* get(const index_t &bi) const
    {
        return valid(bi) ? this->bundle_storage_->get(bi) : nullptr;
    }

    inline size_m_t getSizeM() const
    {
        return size_m_;
    }

    inline size_t getSize() const
    {
        return size_;
    }

    inline size_t getBundleSize() const
    {
        return size_ * 2;
    }

protected:
    const size_t    size_;
    const size_m_t  size_m_;

    virtual inline void updateIndices(const index_t &chunk_index) const override
    {
    }

    virtual inline bool valid(const index_t &index) const override
    {
        auto is_valid = [this,&index](const std::size_t &c) {
            return index[c] >= this->min_bundle_index_[c] && index[c] <= this->max_bundle_index_[c];
        };
        for (std::size_t i=0; i<Dim; ++i)
            if (!is_valid(i))
                return false;
        return true;
    }

    virtual inline void allocateStorage(distribution_storage_t& storage) const override
    {
        storage.template set<cis::option::tags::array_size>(std::forward(this->size_ * 2));
        storage.template set<cis::option::tags::array_offset>(std::forward(this->min_bundle_index_));
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 GenericMap<tags::dynamic_map, Dim, data_t, T, backend_t> : public AbstractMap<tags::dynamic_map, Dim, data_t, T, backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t   = Eigen::aligned_allocator<GenericMap<tags::dynamic_map, Dim, data_t, T, backend_t>>;

    using ConstPtr      = std::shared_ptr<const GenericMap<tags::dynamic_map, Dim, data_t, T, backend_t>>;
    using Ptr           = std::shared_ptr<GenericMap<tags::dynamic_map, Dim, data_t, T, backend_t>>;

    using base_t = AbstractMap<tags::dynamic_map, Dim, data_t, T, backend_t>;
    using typename base_t::pose_t;
    using typename base_t::transform_t;
    using typename base_t::point_t;
    using typename base_t::pointcloud_t;
    using typename base_t::index_t;
    using typename base_t::index_list_t;
    using typename base_t::distribution_t;
    using typename base_t::distribution_storage_t;
    using typename base_t::distribution_storage_ptr_t;
    using typename base_t::distribution_storage_array_t;
    using typename base_t::distribution_bundle_t;
    using typename base_t::distribution_const_bundle_t;
    using typename base_t::distribution_bundle_storage_t;
    using typename base_t::distribution_bundle_storage_ptr_t;

    inline GenericMap(const T resolution) :
        GenericMap(pose_t::identity(), resolution)
    {
    }

    inline GenericMap(const pose_t &origin,
                      const T      &resolution) :
        base_t(origin, resolution,
               utility::create<int,Dim>(std::numeric_limits<int>::max()),
               utility::create<int,Dim>(std::numeric_limits<int>::min()))
    {
    }

    inline GenericMap(const pose_t  &origin,
                      const T       &resolution,
                      const index_t &min_index,
                      const index_t &max_index,
                      const distribution_bundle_storage_ptr_t &bundles,
                      const distribution_storage_array_t      &storage) :
        base_t(origin, resolution, min_index, max_index, bundles, storage)
    {
    }

    inline bool empty() const
    {
        return this->min_bundle_index_[0] == std::numeric_limits<int>::max();
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return this->getAllocate(bi);
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return this->getAllocate(bi);
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        const index_t bi = this->toBundleIndex(p);
        return this->getAllocate(bi);
    }

    inline const distribution_bundle_t* get(const point_t &p) const
    {
        const index_t bi = this->toBundleIndex(p);
        return this->bundle_storage_->get(bi);
    }

    inline const distribution_bundle_t* get(const index_t &bi) const
    {
        return this->bundle_storage_->get(bi);
    }

protected:
    virtual inline void updateIndices(const index_t &chunk_index) const override
    {
        this->min_bundle_index_ = std::min(this->min_bundle_index_, chunk_index);
        this->max_bundle_index_ = std::max(this->max_bundle_index_, chunk_index);
    }

    virtual inline bool valid(const index_t &index) const override
    {
        return true;
    }

    virtual inline void allocateStorage(distribution_storage_t& storage) const override
    {
    }
};
}
}

#endif // CSLIBS_NDT_MAP_GENERIC_MAP_HPP
