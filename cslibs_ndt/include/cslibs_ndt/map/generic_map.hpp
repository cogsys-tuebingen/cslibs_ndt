#ifndef CSLIBS_NDT_MAP_GENERIC_MAP_HPP
#define CSLIBS_NDT_MAP_GENERIC_MAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_ndt/map/traits.hpp>
#include <cslibs_ndt/utility/utility.hpp>

#include <cslibs_math/common/array.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/operations/clustering.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = option_t::default_backend_t>
class EIGEN_ALIGN16 GenericMap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t   = Eigen::aligned_allocator<GenericMap<option_t,Dim,data_t,T,backend_t>>;

    using ConstPtr      = std::shared_ptr<const GenericMap<option_t,Dim,data_t,T,backend_t>>;
    using Ptr           = std::shared_ptr<GenericMap<option_t,Dim,data_t,T,backend_t>>;

    using pose_2d_t     = typename traits<Dim,T>::pose_2d_t;
    using pose_t        = typename traits<Dim,T>::pose_t;
    using transform_t   = typename traits<Dim,T>::transform_t;
    using point_t       = typename traits<Dim,T>::point_t;
    using pointcloud_t  = typename traits<Dim,T>::pointcloud_t;
    using index_t       = std::array<int,Dim>;

    static constexpr std::size_t bin_count  = utility::two_pow(Dim);

    using index_list_t                      = std::array<index_t, bin_count>;
    using distribution_t                    = data_t<T,Dim>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, backend_t>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, bin_count>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, bin_count>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, bin_count>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, backend_t>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;

    using neighborhood_t = cis::operations::clustering::GridNeighborhoodStatic<std::tuple_size<index_t>::value, Dim>;

    template <std::size_t DD>;
    using vector_t = cslibs_math::linear::Vector<T,DD>;

    inline GenericMap(const GenericMap &other);
    inline GenericMap(GenericMap &&other);

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    template <std::size_t... counter>
    inline point_t getMin(utility::integer_sequence<std::size_t,counter...>) const
    {
        auto at = [](const std::size_t &c) {
            return (min_bundle_index_[c] + 1) * bundle_resolution_;
        };
        return point_t(at(counter)...);
    }
    inline point_t getMin() const
    {
        return getMin(utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    template <std::size_t... counter>
    inline point_t getMax(utility::integer_sequence<std::size_t,counter...>) const
    {
        auto at = [](const std::size_t &c) {
            return (max_bundle_index_[c] + 1) * bundle_resolution_;
        };
        return point_t(at(counter)...);
    }
    inline point_t getMax() const
    {
        return getMax(utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
    }

    /**
     * @brief Get the origin.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
        origin.translation() = getMin();
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

    inline index_t getMinBundleIndex() const
    {
        return min_bundle_index_;
    }

    inline index_t getMaxBundleIndex() const
    {
        return max_bundle_index_;
    }

    inline T getBundleResolution() const
    {
        return bundle_resolution_;
    }

    inline T getResolution() const
    {
        return resolution_;
    }

    inline T getHeight() const
    {
        return (max_bundle_index_[1] - min_bundle_index_[1] + 1) * bundle_resolution_;
    }

    inline T getWidth() const
    {
        return (max_bundle_index_[0] - min_bundle_index_[0] + 1) * bundle_resolution_;
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const;
    inline distribution_bundle_t* getDistributionBundle(const index_t &bi);
    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const;
    inline const distribution_bundle_t* get(const point_t &p) const;
    inline const distribution_bundle_t* get(const index_t &bi) const;

    inline distribution_storage_array_t const & getStorages() const
    {
        return storage_;
    }

    template <typename Fn>
    inline void traverse(const Fn& function) const
    {
        return bundle_storage_->traverse(function);
    }

    inline void getBundleIndices(std::vector<index_t> &indices) const
    {
        auto add_index = [&indices](const index_t &i, const distribution_bundle_t &) {
            indices.emplace_back(i);
        };
        bundle_storage_->traverse(add_index);
    }

    inline virtual bool validate(const pose_2d_t &p_w_2d) const
    {
        const point_t p_w = toPoint(p_w_2d.translation());
        return valid(toBundleIndex(p_w));
    }

    inline void allocatePartiallyAllocatedBundles()
    {
        std::vector<index_t> bis;
        getBundleIndices(bis);

        static constexpr neighborhood_t grid{};

        for (const index_t &bi : bis) {
            const distribution_bundle_t *bundle = bundle_storage_->get(bi);
            bool do_expand = expandBundle(bundle);
            if (do_expand) {
                grid.visit([this, &bi](neighborhood_t::offset_t o) {
                    getAllocate(add(bi,o));
                });
            }
        }
    }

    inline std::size_t getByteSize() const
    {
        std::size_t size = sizeof(*this) + bundle_storage_->byte_size();
        for (auto &storage : storage_)
            size += storage->byte_size();
        return size;
    }

protected:
    const T                                    resolution_;
    const T                                    bundle_resolution_;
    const T                                    bundle_resolution_inv_;
    const transform_t                          w_T_m_;
    const transform_t                          m_T_w_;

    mutable index_t                            min_bundle_index_;
    mutable index_t                            max_bundle_index_;
    mutable distribution_storage_array_t       storage_;
    mutable distribution_bundle_storage_ptr_t  bundle_storage_;

    inline static distribution_t* getAllocate(const distribution_storage_ptr_t &s,
                                              const index_t &i) const
    {
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
    {
        auto get_allocate = [this](const index_t &bi) {
            distribution_bundle_t *bundle = bundle_storage_->get(bi);

            auto allocate_bundle = [this, &bi]() {
                static constexpr index_list_t indices = utility::generate_indices<index_list_t,Dim>(bi);

                distribution_bundle_t b;
                std::size_t id = 0;
                for (const auto& index : indices) {
                    b[id] = getAllocate(storage_[id], index);
                    ++id;
                }
                updateIndices(bi);
                return &(bundle_storage_->insert(bi, b));
            };
            return bundle ? bundle : allocate_bundle();
        };

        return get_allocate(bi);
    }

    inline void updateIndices(const index_t &chunk_index) const;
    inline bool valid(const index_t &index) const;

    template <typename std::size_t... counter>
    inline index_t toBundleIndex(const point_t &p_m, utility::integer_sequence<std::size_t,counter...>)
    {
        auto to_bundle_index = [&p_m](const std::size_t &c) {
            return static_cast<int>(std::floor(p_m(c) * bundle_resolution_inv_));
        };
        return {to_bundle_index(counter)...};
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return toBundleIndex(p_m, utility::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
    }

    inline bool toBundleIndex(const point_t &p_w,
                              index_t &index) const
    {
        index = toBundleIndex(p_w);
        return valid(index);
    }

    inline void allocateStorage(distribution_storage_t& storage) const;

protected:
    static inline bool expandDistribution(const distribution_t* d) const;

    template <typename std::size_t... counter>
    static inline bool expandBundle(const distribution_bundle_t *bundle, utility::integer_sequence<std::size_t,counter...>) const
    {
        auto expand_bundle = [&bundle](const std::size_t &c) {
            return expandDistribution(bundle->at(c));
        };
        return utility::merge<utility::operations::bool_or>(expand_bundle(counter)...);
    }

    static inline bool expandBundle(const distribution_bundle_t *bundle) const
    {
        return bundle && expandBundle(bundle, utility::make_integer_sequence<std::size_t,bin_count>{});
    }

    template <typename std::size_t... counter>
    static inline index_t add(const index_t &bi, const neighborhood_t::offset_t &o, utility::integer_sequence<std::size_t,counter...>)
    {
        auto add = [&bi,&o](const std::size_t &c) {
            return bi[c] + o[c];
        };
        return {add(counter)...};
    }

    static inline index_t add(const index_t &bi, const neighborhood_t::offset_t &o) const
    {
        return add(bi, o, utility::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
    }

    template <std::size_t DD, typename std::size_t... counter>
    static inline point_t toPoint(vector_t<DD> p, utility::integer_sequence<std::size_t,counter...>)
    {
        auto at = [&p](const std::size_t &c) {
            return (c >= DD) ? T(0) : p(c);
        };
        return point_t(at(counter)...);
    }

    template <std::size_t DD>
    static inline point_t toPoint(vector_t<DD> p)
    {
        return toPoint(p, utility::make_integer_sequence<std::size_t, point_t::Dimension>{});
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 GenericMap<tags::static_map, Dim, data_t, T, backend_t>
{
public:
    using size_t   = std::array<std::size_t,Dim>;
    using size_m_t = std::array<T,Dim>;

    inline GenericMap(const pose_t  &origin,
                      const T       &resolution,
                      const size_t  &size,
                      const index_t &min_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        min_bundle_index_(min_bundle_index),
        storage_(utility::create<distribution_storage_t,Dim>()),
        bundle_storage_(new distribution_bundle_storage_t)
    {
        for (std::size_t i=0; i<Dim; ++i) {
            size_m_[i] = (size[i] + 1) * resolution;
            max_bundle_index_[i] = min_bundle_index[i] + static_cast<int>(size[i] * 2) -1;
        }

        index_t offset;
        for (std::size_t i=0; i<Dim; ++i)
            offset[i] = cslibs_math::common::div<int>(min_bundle_index[i], 2);
        storage_[0]->template set<cis::option::tags::array_size>(std::forward(size));
        storage_[0]->template set<cis::option::tags::array_offset>(std::forward(offset));

        for(std::size_t i=1 ; i<Dim; ++i) {
            storage_[i]->template set<cis::option::tags::array_size>(std::forward(size+1));
            storage_[i]->template set<cis::option::tags::array_offset>(offset);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(std::forward(size*2));
        bundle_storage_->template set<cis::option::tags::array_offset>(std::forward(min_bundle_index));
    }

    inline GenericMap(const pose_t &origin,
                      const T      &resolution,
                      const size_t &size,
                      const distribution_bundle_storage_ptr_t &bundles,
                      const distribution_storage_array_t      &storage,
                      const index_t                           &min_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        min_bundle_index_(min_bundle_index),
        storage_(storage),
        bundle_storage_(bundles)
    {
        for (std::size_t i=0; i<Dim; ++i) {
            size_m_[i] = (size[i] + 1) * resolution;
            max_bundle_index_[i] = min_bundle_index[i] + static_cast<int>(size[i] * 2) -1;
        }
    }

    inline GenericMap(const GenericMap &other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(other.w_T_m_),
        m_T_w_(other.m_T_w_),
        size_(other.size_),
        size_m_(other.size_m_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(utility::create<distribution_storage_t,Dim>(other.storage_)),
        bundle_storage_(new distribution_bundle_storage_t(*other.bundle_storage_))
    {
    }

    inline GenericMap(GenericMap &&other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(std::move(other.w_T_m_)),
        m_T_w_(std::move(other.m_T_w_)),
        size_(other.size_),
        size_m_(other.size_m_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(other.storage_),
        bundle_storage_(other.bundle_storage_)
    {
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return valid(bi) ? getAllocate(bi) : nullptr;
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return valid(bi) ? getAllocate(bi) : nullptr;
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        index_t bi;
        if(!toBundleIndex(p, bi))
            return nullptr;

        return getAllocate(bi);
    }

    inline const distribution_bundle_t* get(const point_t &p) const
    {
        index_t bi;
        if(!toBundleIndex(p, bi))
            return nullptr;

        return bundle_storage_->get(bi);
    }

    inline const distribution_bundle_t* get(const index_t &bi) const
    {
        return valid(bi) ? bundle_storage_->get(bi) : nullptr;
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

    inline void updateIndices(const index_t &chunk_index) const
    {
    }

    inline bool valid(const index_t &index) const
    {
        auto is_valid = [&index](const std::size_t &c) {
            return index[c] >= min_bundle_index_[c] && index[c] <= max_bundle_index_[c];
        };
        for (std::size_t i=0; i<Dim; ++i)
            if (!is_valid(i))
                return false;
        return true;
    }

    inline void allocateStorage(distribution_storage_t& storage) const
    {
        storage.template set<cis::option::tags::array_size>(std::forward(size_ * 2));
        storage.template set<cis::option::tags::array_offset>(std::forward(min_bundle_index_));
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 GenericMap<tags::dynamic_map, Dim, data_t, T, backend_t>
{
public:
    inline GenericMap(const T resolution) :
        GenericMap(pose_t::identity(),
                   resolution)
    {
    }

    inline GenericMap(const pose_t &origin,
                      const T      &resolution) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_bundle_index_(utility::create<int,Dim>(std::numeric_limits<int>::max())),
        max_bundle_index_(utility::create<int,Dim>(std::numeric_limits<int>::min())),
        storage_(utility::create<distribution_storage_t,Dim>()),
        bundle_storage_(new distribution_bundle_storage_t)
    {
    }

    inline GenericMap(const pose_t  &origin,
                      const T       &resolution,
                      const index_t &min_index,
                      const index_t &max_index,
                      const distribution_bundle_storage_ptr_t &bundles,
                      const distribution_storage_array_t      &storage) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_bundle_index_(min_index),
        max_bundle_index_(max_index),
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline GenericMap(const GenericMap &other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(other.w_T_m_),
        m_T_w_(other.m_T_w_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(utility::create<distribution_storage_t,Dim>(other.storage_)),
        bundle_storage_(new distribution_bundle_storage_t(*other.bundle_storage_))
    {
    }

    inline GenericMap(GenericMap &&other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(std::move(other.w_T_m_)),
        m_T_w_(std::move(other.m_T_w_)),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(other.storage_),
        bundle_storage_(other.bundle_storage_)
    {
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return getAllocate(bi);
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return getAllocate(bi);
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        return getAllocate(bi);
    }

    inline const distribution_bundle_t* get(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        return bundle_storage_->get(bi);
    }

    inline const distribution_bundle_t* get(const index_t &bi) const
    {
        return bundle_storage_->get(bi);
    }

protected:
    inline void updateIndices(const index_t &chunk_index) const
    {
        min_bundle_index_ = std::min(min_bundle_index_, chunk_index);
        max_bundle_index_ = std::max(max_bundle_index_, chunk_index);
    }

    inline bool valid(const index_t &index) const
    {
        return true;
    }

    inline void allocateStorage(distribution_storage_t& storage) const
    {
    }
};
}
}

#endif // CSLIBS_NDT_MAP_GENERIC_MAP_HPP
