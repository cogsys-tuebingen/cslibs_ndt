#ifndef CSLIBS_NDT_MAP_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_ndt/map/utility.hpp>
#include <cslibs_indexed_storage/storage.hpp>
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

    inline distribution_t* getAllocate(const distribution_storage_ptr_t &s,
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

    inline int toBundleIndex(const point_t &p_m, const std::size_t &counter)
    {
        return static_cast<int>(std::floor(p_m(counter) * bundle_resolution_inv_));
    }

    template <typename std::size_t... counter>
    inline index_t toBundleIndex(const point_t &p_m, detail::integer_sequence<std::size_t,counter...>)
    {
        return {toBundleIndex(p_m,counter)...};
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return toBundleIndex(p_m, detail::make_integer_sequence<std::size_t,std::tuple_size<index_t>::value>{});
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

    inline bool valid(const index_t &index, const std::size_t &counter)
    {
        return index[counter] >= min_bundle_index_[counter] && index[counter] <= max_bundle_index_[counter];
    }

    template <std::size_t... counter>
    inline bool valid(const index_t &index, detail::integer_sequence<std::size_t,counter...> counts)
    {
        return detail::merge_and(valid(index, counter)...);
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
