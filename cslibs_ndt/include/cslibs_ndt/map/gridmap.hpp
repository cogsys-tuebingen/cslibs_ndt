#ifndef CSLIBS_NDT_MAP_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_ndt/map/utility.hpp>
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
class EIGEN_ALIGN16 Gridmap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t   = Eigen::aligned_allocator<Gridmap<option_t,Dim,data_t,T,backend_t>>;

    using ConstPtr      = std::shared_ptr<const Gridmap<option_t,Dim,data_t,T,backend_t>>;
    using Ptr           = std::shared_ptr<Gridmap<option_t,Dim,data_t,T,backend_t>>;

    using pose_t        = typename traits<Dim,T>::pose_t;
    using transform_t   = typename traits<Dim,T>::transform_t;
    using point_t       = typename traits<Dim,T>::point_t;
    using pointcloud_t  = typename traits<Dim,T>::pointcloud_t;
    using index_t       = std::array<int,Dim>;

    static constexpr std::size_t bin_count  = detail::two_pow(Dim);

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
                static constexpr index_list_t indices = detail::generate_indices<index_list_t,Dim>(bi);

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

    inline void updateIndices(const index_t &chunk_index) const
    {
        // default: do nothing,
        // dynamic maps have to update min and max index
    }

    template <typename std::size_t... counter>
    inline index_t toBundleIndex(const point_t &p_m, detail::integer_sequence<std::size_t,counter...>)
    {
        auto to_bundle_index = [&p_m](const std::size_t &c) {
            return static_cast<int>(std::floor(p_m(c) * bundle_resolution_inv_));
        };
        return {to_bundle_index(counter)...};
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return toBundleIndex(p_m, detail::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
    }

protected:
    static inline bool expandDistribution(const distribution_t* d) const;

    template <typename std::size_t... counter>
    static inline bool expandBundle(const distribution_bundle_t *bundle, detail::integer_sequence<std::size_t,counter...>) const
    {
        auto expand_bundle = [&bundle](const std::size_t &c) {
            return expandDistribution(bundle->at(c));
        };
        return detail::merge<detail::bool_or>(expand_bundle(counter)...);
    }

    static inline bool expandBundle(const distribution_bundle_t *bundle) const
    {
        return bundle && expandBundle(bundle, detail::make_integer_sequence<std::size_t,bin_count>{});
    }

    template <typename std::size_t... counter>
    static inline index_t add(const index_t &bi, const neighborhood_t::offset_t &o, detail::integer_sequence<std::size_t,counter...>)
    {
        auto add = [&bi,&o](const std::size_t &c) {
            return bi[c] + o[c];
        };
        return {add(counter)...};
    }

    static inline index_t add(const index_t &bi, const neighborhood_t::offset_t &o) const
    {
        return add(bi, o, detail::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 Gridmap<tags::static_map, Dim, data_t, T, backend_t>
{
public:
    using size_t   = std::array<std::size_t,Dim>;
    using size_m_t = std::array<T,Dim>;

protected:
    const size_t    size_;
    const size_m_t  size_m_;

    template <std::size_t... counter>
    inline bool valid(const index_t &index, detail::integer_sequence<std::size_t,counter...> counts)
    {
        auto valid = [&index](const std::size_t &c) {
            return index[c] >= min_bundle_index_[c] && index[c] <= max_bundle_index_[c];
        };
        return detail::merge<detail::bool_and>(valid(counter)...);
    }

    inline bool valid(const index_t &index) const
    {
        return valid(index,detail::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
    }

    inline bool toBundleIndex(const point_t &p_w,
                              index_t &index) const
    {
        index = toBundleIndex(p_w);
        return valid(index);
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> data_t,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 Gridmap<tags::dynamic_map, Dim, data_t, T, backend_t>
{
protected:
    inline void updateIndices(const index_t &chunk_index) const
    {
        min_bundle_index_ = std::min(min_bundle_index_, chunk_index);
        max_bundle_index_ = std::max(max_bundle_index_, chunk_index);
    }
};
}
}

#endif // CSLIBS_NDT_MAP_GRIDMAP_HPP
