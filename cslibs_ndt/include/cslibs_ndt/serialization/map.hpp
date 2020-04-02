#ifndef CSLIBS_NDT_SERIALIZATION_MAP_HPP
#define CSLIBS_NDT_SERIALIZATION_MAP_HPP

#include <cslibs_ndt/map/map.hpp>
#include <cslibs_ndt/serialization/loader.hpp>

#include <cslibs_ndt/serialization/filesystem.hpp>
#include <cslibs_ndt/serialization/storage.hpp>

#include <cslibs_math_2d/serialization/transform.hpp>
#include <cslibs_math_3d/serialization/transform.hpp>
#include <cslibs_math/serialization/array.hpp>

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <thread>
#include <atomic>

namespace cslibs_ndt {
namespace serialization {

template <map::tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = map::tags::default_types<option_t>::template default_backend_t>
inline bool saveBinary(const cslibs_ndt::map::Map<option_t,Dim,data_t,T,backend_t> &map,
                       const std::string &path)
{
    using path_t     = boost::filesystem::path;
    using map_t      = cslibs_ndt::map::Map<option_t,Dim,data_t,T,backend_t>;
    using paths_t    = std::array<path_t, map_t::bin_count>;
    using binary_t   = cslibs_ndt::binary<data_t, T, Dim, Dim, backend_t>;
    using storages_t = typename map_t::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::create_directory(path_root))
        return false;

    /// step two: identify subfolders
    paths_t paths;
    for (std::size_t i=0; i<map_t::bin_count; ++i)
        paths[i] = path_root / path_t("store_" + std::to_string(i) + ".bin");

    /// step three: we have our filesystem, now we write out the distributions file by file
    /// meta file
    const path_t path_file = path_root / path_t("map.bin");
    if (!header<option_t,Dim,data_t,T,backend_t>::save(map,path_file))
        return false;

    /// step four: write out the storages
    storages_t storages = map.getStorages();

    std::array<std::thread, map_t::bin_count> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < map_t::bin_count; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::save(storages[i], paths[i]);
        });
    for (std::size_t i = 0 ; i < map_t::bin_count; ++i)
        if (threads[i].joinable())
            threads[i].join();

    return success;
}

template <map::tags::option option_t,
          std::size_t Dim,
          template <typename,std::size_t> class data_t,
          typename T,
          template <typename, typename, typename...> class backend_t = map::tags::default_types<option_t>::template default_backend_t>
inline bool loadBinary(const std::string &path,
                       typename cslibs_ndt::map::Map<option_t,Dim,data_t,T,backend_t>::Ptr &map)
{
    using path_t            = boost::filesystem::path;
    using map_t             = cslibs_ndt::map::Map<option_t,Dim,data_t,T,backend_t>;
    using paths_t           = std::array<path_t, map_t::bin_count>;
    using bundle_storage_t  = typename map_t::distribution_bundle_storage_t;
    using storages_t        = typename map_t::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::check_directory(path_root))
        return false;

    /// step two: identify subfolders
    paths_t paths;
    for (std::size_t i = 0 ; i < map_t::bin_count ; ++i)
        paths[i] = path_root / path_t("store_" + std::to_string(i) + ".bin");

    /// step three: we have our filesystem, now we can load distributions file by file
    for (std::size_t i = 0 ; i < map_t::bin_count ; ++i)
        if (!cslibs_ndt::common::serialization::check_file(paths[i]))
            return false;

    std::shared_ptr<bundle_storage_t> bundles(new bundle_storage_t);
    storages_t storages;

    /// load meta data
    path_t path_file = path_root / path_t("map.yaml");
    loader<option_t,Dim,data_t,T,backend_t>* l;

    if (cslibs_ndt::common::serialization::check_file(path_file)) {
        YAML::Node n = YAML::LoadFile((path_root / path_file).string());
        l = header<option_t,Dim,data_t,T,backend_t>::load(n);
    } else {
        path_file = path_root / path_t("map.bin");
        if (!cslibs_ndt::common::serialization::check_file(path_file))
            return false;
        l = header<option_t,Dim,data_t,T,backend_t>::load(path_file);
        if (!l)
            return false;
    }

    std::array<std::thread, map_t::bin_count> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < map_t::bin_count ; ++i) {
        const path_t path = paths[i];
        threads[i] = std::thread([&l, &storages, path, i, &success](){
            success = success && l->load(i, path, storages[i]);
        });
    }
    for (std::size_t i = 0 ; i < map_t::bin_count ; ++i)
        threads[i].join();

    if (!success)
        return false;

    l->allocateBundles(bundles, storages);
    l->createMap(bundles, storages, map);

    delete l;
    return true;
}

}
}

#endif // CSLIBS_NDT_SERIALIZATION_MAP_HPP
