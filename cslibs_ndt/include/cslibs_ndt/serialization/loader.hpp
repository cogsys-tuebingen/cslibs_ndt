#ifndef CSLIBS_NDT_SERIALIZATION_LOADER_HPP
#define CSLIBS_NDT_SERIALIZATION_LOADER_HPP

#include <cslibs_ndt/map/map.hpp>

#include <cslibs_ndt/serialization/filesystem.hpp>
#include <cslibs_ndt/serialization/storage.hpp>

#include <thread>
#include <atomic>

namespace cslibs_ndt {
namespace serialization {

template <map::tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct loader {};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct loader<cslibs_ndt::map::tags::static_map,Dim,data_t,T,backend_t,dynamic_backend_t> {
    using map_t             = cslibs_ndt::map::Map<cslibs_ndt::map::tags::static_map,Dim,data_t,T,backend_t,dynamic_backend_t>;
    using index_t           = typename map_t::index_t;
    using pose_t            = typename map_t::pose_t;
    using size_t            = typename map_t::size_t;
    using bundle_storage_t  = typename map_t::distribution_bundle_storage_t;
    using storage_t         = typename map_t::distribution_storage_ptr_t;
    using storages_t        = typename map_t::distribution_storage_array_t;
    using binary_t          = cslibs_ndt::binary<data_t, T, Dim, Dim, backend_t>;
    using path_t            = boost::filesystem::path;
    using paths_t           = std::array<path_t, map_t::bin_count>;

    inline loader(const pose_t &pose,
                  const T& resolution,
                  const size_t &size,
                  const index_t &min_index,
                  const std::vector<index_t> &indices) :
        pose_(pose), resolution_(resolution), size_(size), min_index_(min_index), indices_(indices)
    {
    }

    inline bool load(const std::size_t i, const path_t path, storage_t &storage) const
    {
        const std::size_t off = (i > 1ul) ? 1ul : 0ul;
        index_t offset;
        for (std::size_t i=0; i<Dim; ++i)
            offset[i] = cslibs_math::common::div<int>(min_index_[i], 2);
        return binary_t::load(path, storage, size_ + off, offset);
    }

    inline void allocateBundles(const std::shared_ptr<bundle_storage_t>& bundles,
                                const storages_t& storages) const
    {
        bundles->template set<cslibs_indexed_storage::option::tags::array_size>(size_ * 2ul);
        bundles->template set<cslibs_indexed_storage::option::tags::array_offset>(min_index_);

        auto allocate_bundle = [&storages, &bundles](const index_t &bi) {
            typename map_t::distribution_bundle_t b;
            const typename map_t::index_list_t indices =
                    utility::generate_indices<typename map_t::index_list_t,Dim>(bi);
            for (std::size_t i = 0 ; i < map_t::bin_count ; ++i)
                b[i] = storages[i]->get(indices[i]);
            bundles->insert(bi, b);
        };
        for (const index_t &index : indices_)
            allocate_bundle(index);
    }

    inline void setMap(typename map_t::Ptr &map,
                       const std::shared_ptr<bundle_storage_t>& bundles,
                       const storages_t& storages) const {
        map.reset(new map_t(pose_,
                            resolution_,
                            size_,
                            bundles,
                            storages,
                            min_index_));
    }

    const pose_t pose_;
    const T resolution_;
    const size_t size_;
    const index_t min_index_;
    const std::vector<index_t> indices_;
};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct loader<cslibs_ndt::map::tags::dynamic_map,Dim,data_t,T,backend_t,dynamic_backend_t> {
    using map_t             = cslibs_ndt::map::Map<cslibs_ndt::map::tags::dynamic_map,Dim,data_t,T,backend_t,dynamic_backend_t>;
    using index_t           = typename map_t::index_t;
    using pose_t            = typename map_t::pose_t;
    using bundle_storage_t  = typename map_t::distribution_bundle_storage_t;
    using storage_t         = typename map_t::distribution_storage_ptr_t;
    using storages_t        = typename map_t::distribution_storage_array_t;
    using binary_t          = cslibs_ndt::binary<data_t, T, Dim, Dim, backend_t>;
    using path_t            = boost::filesystem::path;
    using paths_t           = std::array<path_t, map_t::bin_count>;

    inline loader(const pose_t &pose,
                  const T& resolution,
                  const index_t &min_index,
                  const index_t &max_index,
                  const std::vector<index_t> &indices) :
        pose_(pose), resolution_(resolution), min_index_(min_index), max_index_(max_index), indices_(indices)
    {
    }

    inline bool load(const std::size_t i, const path_t path, storage_t &storage) const
    {
        return binary_t::load(path, storage);
    }

    inline void allocateBundles(const std::shared_ptr<bundle_storage_t>& bundles,
                                const storages_t& storages) const
    {
        auto allocate_bundle = [&storages, &bundles](const index_t &bi) {
            typename map_t::distribution_bundle_t b;
            const typename map_t::index_list_t indices =
                    utility::generate_indices<typename map_t::index_list_t,Dim>(bi);
            for (std::size_t i = 0 ; i < map_t::bin_count ; ++i)
                b[i] = storages[i]->get(indices[i]);
            bundles->insert(bi, b);
        };
        for (const index_t &index : indices_)
            allocate_bundle(index);
    }

    inline void setMap(typename map_t::Ptr &map,
                       const std::shared_ptr<bundle_storage_t>& bundles,
                       const storages_t& storages) const {
        map.reset(new map_t(pose_,
                            resolution_,
                            min_index_,
                            max_index_,
                            bundles,
                            storages));
    }

    const pose_t pose_;
    const T resolution_;
    const index_t min_index_;
    const index_t max_index_;
    const std::vector<index_t> indices_;
};


template <map::tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct header {};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct header<cslibs_ndt::map::tags::static_map,Dim,data_t,T,backend_t,dynamic_backend_t> {
    using map_t    = cslibs_ndt::map::Map<cslibs_ndt::map::tags::static_map,Dim,data_t,T,backend_t,dynamic_backend_t>;
    using index_t  = typename map_t::index_t;
    using pose_t   = typename map_t::pose_t;
    using size_t   = typename map_t::size_t;
    using loader_t = loader<cslibs_ndt::map::tags::static_map,Dim,data_t,T,backend_t,dynamic_backend_t>;

    static inline void write(const typename map_t::Ptr &map, YAML::Node &n)
    {
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"]     = map->getInitialOrigin();
        n["resolution"] = map->getResolution();
        n["size"]       = map->getSize();
        n["min_index"]  = map->getMinBundleIndex();
        n["bundles"]    = indices;
    }

    static inline loader_t load(const YAML::Node& n) {
        const pose_t               origin     = n["origin"].as<pose_t>();
        const T                    resolution = n["resolution"].as<T>();
        const size_t               size       = n["size"].as<size_t>();
        const index_t              min_index  = n["min_index"].as<index_t>();
        const std::vector<index_t> indices    = n["bundles"].as<std::vector<index_t>>();

        return loader_t(origin,resolution,size,min_index,indices);
    }
};

template <std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
struct header<cslibs_ndt::map::tags::dynamic_map,Dim,data_t,T,backend_t,dynamic_backend_t> {
    using map_t    = cslibs_ndt::map::Map<cslibs_ndt::map::tags::dynamic_map,Dim,data_t,T,backend_t,dynamic_backend_t>;
    using index_t  = typename map_t::index_t;
    using pose_t   = typename map_t::pose_t;
    using loader_t = loader<cslibs_ndt::map::tags::dynamic_map,Dim,data_t,T,backend_t,dynamic_backend_t>;

    static inline void write(const typename map_t::Ptr &map, YAML::Node &n)
    {
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"]     = map->getInitialOrigin();
        n["resolution"] = map->getResolution();
        n["min_index"]  = map->getMinBundleIndex();
        n["max_index"]  = map->getMaxBundleIndex();
        n["bundles"]    = indices;
    }

    static inline loader_t load(const YAML::Node& n) {
        const pose_t               origin     = n["origin"].as<pose_t>();
        const T                    resolution = n["resolution"].as<T>();
        const index_t              min_index  = n["min_index"].as<index_t>();
        const index_t              max_index  = n["max_index"].as<index_t>();
        const std::vector<index_t> indices    = n["bundles"].as<std::vector<index_t>>();

        return loader_t(origin,resolution,min_index,max_index,indices);
    }
};

}
}

#endif // CSLIBS_NDT_SERIALIZATION_LOADER_HPP
