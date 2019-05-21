#ifndef CSLIBS_NDT_MAP_ABSTRACT_MAP_HPP
#define CSLIBS_NDT_MAP_ABSTRACT_MAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_ndt/map/traits.hpp>
#include <cslibs_ndt/common/bundle.hpp>
#include <cslibs_ndt/utility/utility.hpp>

#include <cslibs_math/common/array.hpp>
#include <cslibs_math/utility/traits.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/operations/clustering.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = tags::default_types<option_t>::template default_backend_t,
          template <typename, typename, typename...> class dynamic_backend_t = tags::default_types<option_t>::template default_dynamic_backend_t>
class EIGEN_ALIGN16 AbstractMap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t   = Eigen::aligned_allocator<AbstractMap<option_t,Dim,data_t,T,backend_t,dynamic_backend_t>>;

    using ConstPtr      = std::shared_ptr<const AbstractMap<option_t,Dim,data_t,T,backend_t,dynamic_backend_t>>;
    using Ptr           = std::shared_ptr<AbstractMap<option_t,Dim,data_t,T,backend_t,dynamic_backend_t>>;

    using pose_2d_t     = typename traits<Dim,T>::pose_2d_t;
    using pose_t        = typename traits<Dim,T>::pose_t;
    using transform_t   = typename traits<Dim,T>::transform_t;
    using point_t       = typename traits<Dim,T>::point_t;
    using pointcloud_t  = typename traits<Dim,T>::pointcloud_t;
    using index_t       = std::array<int,Dim>;

    static constexpr std::size_t bin_count  = utility::two_pow(Dim);
    static constexpr T div_count = cslibs_math::utility::traits<T>::One / static_cast<T>(bin_count);

    using index_list_t                      = std::array<index_t, bin_count>;
    using distribution_t                    = data_t<T,Dim>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, backend_t>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, bin_count>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, bin_count>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, bin_count>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, backend_t>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;
    using dynamic_distribution_storage_t    = cis::Storage<distribution_t, index_t, dynamic_backend_t>;

    using neighborhood_t = cis::operations::clustering::GridNeighborhoodStatic<std::tuple_size<index_t>::value, 3>;

    template <std::size_t DD>
    using vector_t = cslibs_math::linear::Vector<T,DD>;

    inline AbstractMap(const pose_t  &origin,
                       const T       &resolution,
                       const index_t &min_bundle_index,
                       const index_t &max_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(cslibs_math::utility::traits<T>::Half * resolution_),
        bundle_resolution_inv_(cslibs_math::utility::traits<T>::One / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_bundle_index_(min_bundle_index),
        max_bundle_index_(max_bundle_index),
        storage_(utility::create<distribution_storage_t,bin_count>()),
        bundle_storage_(new distribution_bundle_storage_t)
    {
    }

    inline AbstractMap(const pose_t  &origin,
                       const T       &resolution,
                       const index_t &min_bundle_index,
                       const index_t &max_bundle_index,
                       const distribution_bundle_storage_ptr_t &bundles,
                       const distribution_storage_array_t      &storage) :
        resolution_(resolution),
        bundle_resolution_(cslibs_math::utility::traits<T>::Half * resolution_),
        bundle_resolution_inv_(cslibs_math::utility::traits<T>::One / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_bundle_index_(min_bundle_index),
        max_bundle_index_(max_bundle_index),
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline AbstractMap(const AbstractMap &other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(other.w_T_m_),
        m_T_w_(other.m_T_w_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(utility::create<distribution_storage_t,bin_count>(other.storage_)),
        bundle_storage_(new distribution_bundle_storage_t(*other.bundle_storage_))
    {
    }

    inline AbstractMap(AbstractMap &&other) :
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

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    inline point_t getMin() const
    {
        point_t min;
        for (std::size_t i=0; i<Dim; ++i)
            min(i) = static_cast<T>(min_bundle_index_[i]) * bundle_resolution_;
        return min;
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    inline point_t getMax() const
    {
        point_t max;
        for (std::size_t i=0; i<Dim; ++i)
            max(i) = static_cast<T>(max_bundle_index_[i] + 1) * bundle_resolution_;
        return max;
    }

    /**
     * @brief Get the origin.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
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

    inline void allocatePartiallyAllocatedBundles() const
    {
        std::vector<index_t> bis;
        getBundleIndices(bis);

        static constexpr neighborhood_t grid{};

        for (const index_t &bi : bis) {
            const distribution_bundle_t *bundle = bundle_storage_->get(bi);
            bool do_expand = expandBundle(bundle);
            if (do_expand) {
                grid.visit([this, &bi](typename neighborhood_t::offset_t o) {
                    index_t ii;
                    for (std::size_t i=0; i<Dim; ++i)
                        ii[i] = bi[i] + o[i];
                    if (valid(ii))
                        getAllocate(ii);
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
                                              const index_t &i)
    {
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
    {
        auto get_allocate = [this](const index_t &bi) {
            distribution_bundle_t *bundle = bundle_storage_->get(bi);

            auto allocate_bundle = [this, &bi]() {
                const index_list_t indices = utility::generate_indices<index_list_t,Dim>(bi);

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

    virtual void updateIndices(const index_t &chunk_index) const = 0;
    virtual bool valid(const index_t &index) const = 0;

    inline index_t toBundleIndex(const point_t &p_w,
                                 point_t &p_m) const
    {
        p_m = m_T_w_ * p_w;
        index_t retval;
        for (std::size_t i=0; i<Dim; ++i)
            retval[i] = static_cast<int>(std::floor(p_m(i) * bundle_resolution_inv_));
        return retval;
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        point_t p_m;
        return toBundleIndex(p_w, p_m);
    }

    inline bool toBundleIndex(const point_t &p_w,
                              index_t &index) const
    {
        index = toBundleIndex(p_w);
        return valid(index);
    }

    inline bool toBundleIndex(const point_t &p_w,
                              point_t &p_m,
                              index_t &index) const
    {
        index = toBundleIndex(p_w, p_m);
        return valid(index);
    }

    virtual bool expandDistribution(const distribution_t* d) const = 0;

    inline bool expandBundle(const distribution_bundle_t *bundle) const
    {
        if (!bundle)
            return false;
        for (std::size_t i=0; i<bin_count; ++i)
            if (expandDistribution(bundle->at(i)))
                return true;
        return false;
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
}
}

#endif // CSLIBS_NDT_MAP_ABSTRACT_MAP_HPP
