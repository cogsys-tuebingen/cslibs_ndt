#ifndef CSLIBS_NDT_CONVERSION_MAP_HPP
#define CSLIBS_NDT_CONVERSION_MAP_HPP

#include <cslibs_ndt/map/map.hpp>
#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/occupancy_distribution.hpp>

namespace cslibs_ndt {
namespace conversion {

namespace impl {
template <template <typename,std::size_t> class data_t, typename T, std::size_t Dim>
struct convert {};

template <typename T, std::size_t Dim>
struct convert<Distribution,T,Dim> {
    static inline void from(const Distribution<T,Dim>* const& f, Distribution<T,Dim>* const& t)
    {
        t->data() = f->data();
    }
};

template <typename T, std::size_t Dim>
struct convert<OccupancyDistribution,T,Dim> {
    static inline void from(const OccupancyDistribution<T,Dim>* const& f, OccupancyDistribution<T,Dim>* const& t)
    {
        if (f && (f->numFree() > 0 || f->numOccupied() > 0))
            *t = *f;
    }
};
}

template <map::tags::option option_to_t,
          map::tags::option option_from_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_to_t = map::tags::default_types<option_to_t>::template default_backend_t,
          template <typename, typename, typename...> class backend_from_t = map::tags::default_types<option_from_t>::template default_backend_t>
struct convert {};

template <map::tags::option option_from_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_to_t,
          template <typename, typename, typename...> class backend_from_t>
struct convert<map::tags::dynamic_map,option_from_t,Dim,data_t,T,backend_to_t,backend_from_t> {
    using src_map_t = map::Map<option_from_t,Dim,data_t,T,backend_from_t>;
    using dst_map_t = map::Map<map::tags::dynamic_map,Dim,data_t,T,backend_to_t>;

    static inline typename dst_map_t::Ptr from(const typename src_map_t::Ptr& src)
    {
        if (!src)
            return nullptr;

        typename dst_map_t::Ptr dst(new dst_map_t(src->getInitialOrigin(),
                                                  src->getResolution()));

        static constexpr std::size_t bin_count  = utility::two_pow(Dim);
        using index_t = typename src_map_t::index_t;
        using bundle_t = cslibs_ndt::Bundle<data_t<T,Dim>*, bin_count>;
        src->traverse([&dst](const index_t &bi, const bundle_t &b) {
            if (const bundle_t* b_dst = dst->getDistributionBundle(bi)) {
                for (std::size_t i = 0 ; i < bin_count ; ++i)
                    impl::convert<data_t,T,Dim>::from(b.at(i), b_dst->at(i));
            }
        });

        return dst;
    }
};

template <map::tags::option option_from_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_to_t,
          template <typename, typename, typename...> class backend_from_t>
struct convert<map::tags::static_map,option_from_t,Dim,data_t,T,backend_to_t,backend_from_t> {
    using src_map_t = map::Map<option_from_t,Dim,data_t,T,backend_from_t>;
    using dst_map_t = map::Map<map::tags::static_map,Dim,data_t,T,backend_to_t>;

    static inline typename dst_map_t::Ptr from(const typename src_map_t::Ptr& src)
    {
        if (!src)
            return nullptr;

        using index_t = typename src_map_t::index_t;
        const index_t min_distribution_index =
                cslibs_math::common::cast<int>(std::floor(cslibs_math::common::cast<T>(src->getMinBundleIndex()) / 2.0) * 2.0);
        const index_t max_distribution_index =
                cslibs_math::common::cast<int>( std::ceil(cslibs_math::common::cast<T>(src->getMaxBundleIndex()) / 2.0) * 2.0 + 1.0);

        const typename dst_map_t::size_t size =
                cslibs_math::common::cast<std::size_t>(std::ceil(cslibs_math::common::cast<T>(max_distribution_index - min_distribution_index) / 2.0));

        typename dst_map_t::Ptr dst(new dst_map_t(src->getInitialOrigin(),
                                                  src->getResolution(),
                                                  size,
                                                  min_distribution_index));

        static constexpr std::size_t bin_count  = utility::two_pow(Dim);
        using bundle_t = cslibs_ndt::Bundle<data_t<T,Dim>*, bin_count>;
        src->traverse([&dst](const index_t &bi, const bundle_t &b) {
            if (const bundle_t* b_dst = dst->getDistributionBundle(bi)) {
                for (std::size_t i = 0 ; i < bin_count ; ++i)
                    impl::convert<data_t,T,Dim>::from(b.at(i), b_dst->at(i));
            }
        });

        return dst;
    }
};

}
}

#endif // CSLIBS_NDT_CONVERSION_MAP_HPP
