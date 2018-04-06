#ifndef CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt_2d/static_maps/gridmap.hpp>

#include <cslibs_ndt/serialization/filesystem.hpp>
#include <cslibs_ndt/serialization/storage.hpp>

#include <cslibs_math_2d/serialization/transform.hpp>
#include <cslibs_math/serialization/array.hpp>

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <thread>
#include <atomic>

namespace cslibs_ndt_2d {
namespace static_maps {
inline bool saveBinary(const cslibs_ndt_2d::static_maps::Gridmap::Ptr &map,
                       const std::string &path)
{
    using path_t     = boost::filesystem::path;
    using paths_t    = std::array<path_t, 4>;
    using index_t    = cslibs_ndt_2d::static_maps::Gridmap::index_t;
    using storages_t = cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_array_t;
    using binary_t   = cslibs_ndt::binary<cslibs_ndt::Distribution, 2, 2>;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::create_directory(path_root))
        return false;

    /// step two: check if the sub folders can be created
    const paths_t paths = {{path_root / path_t("store_0.bin"),
                            path_root / path_t("store_1.bin"),
                            path_root / path_t("store_2.bin"),
                            path_root / path_t("store_3.bin")}};

    /// step three: we have our filesystem, now we write out the distributions file by file
    /// meta file
    const path_t path_file = path_t("map.yaml");
    {
        std::ofstream out((path_root / path_file).string(), std::fstream::trunc);
        YAML::Emitter yaml(out);
        YAML::Node n;
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"]     = map->getOrigin();
        n["resolution"] = map->getResolution();
        n["size"]       = map->getSize();
        n["bundles"]    = indices;
        yaml << n;
    }

    /// step four: write out the storages
    const storages_t storages = {{map->getStorages()[0],
                                  map->getStorages()[1],
                                  map->getStorages()[2],
                                  map->getStorages()[3]}};

    std::array<std::thread, 4> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < 4 ; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::save(storages[i], paths[i]);
        });
    for (std::size_t i = 0 ; i < 4 ; ++i)
        threads[i].join();

    return success;
}

inline bool loadBinary(const std::string &path,
                       cslibs_ndt_2d::static_maps::Gridmap::Ptr &map)
{
    using path_t           = boost::filesystem::path;
    using paths_t          = std::array<path_t, 4>;
    using index_t          = cslibs_ndt_2d::static_maps::Gridmap::index_t;
    using size_t           = cslibs_ndt_2d::static_maps::Gridmap::size_t;
    using binary_t         = cslibs_ndt::binary<cslibs_ndt::Distribution, 2, 2>;
    using bundle_storage_t = cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_storage_t;
    using storages_t       = cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::check_directory(path_root))
        return false;

    /// step two: identity subfolders
    const paths_t paths = {{path_root / path_t("store_0.bin"),
                            path_root / path_t("store_1.bin"),
                            path_root / path_t("store_2.bin"),
                            path_root / path_t("store_3.bin")}};

    /// step three: we have our filesystem, now we can load distributions file by file
    for (std::size_t i = 0 ; i < 4 ; ++i)
        if (!cslibs_ndt::common::serialization::check_file(paths[i]))
            return false;

    /// load meta data
    path_t  path_file = path_t("map.yaml");

    std::shared_ptr<bundle_storage_t> bundles(new bundle_storage_t);
    storages_t storages;

    YAML::Node n = YAML::LoadFile((path_root / path_file).string());
    const cslibs_math_2d::Transform2d origin     = n["origin"].as<cslibs_math_2d::Transform2d>();
    const double                      resolution = n["resolution"].as<double>();
    const size_t                      size       = n["size"].as<size_t>();
    const std::vector<index_t>        indices    = n["bundles"].as<std::vector<index_t>>();
    bundles->template set<cslibs_indexed_storage::option::tags::array_size>(size[0] * 2, size[1] * 2);

    std::array<std::thread, 4> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < 4 ; ++ i) {
        const int off   = (i > 1) ? 1 : 0;
        const size_t sz = {{size[0] + off, size[1] + off}};
        threads[i] = std::thread([&storages, &paths, i, &sz, &success](){
            success = success && binary_t::load(paths[i], storages[i], sz);
        });
    }
    for (std::size_t i = 0 ; i < 4 ; ++ i)
        threads[i].join();

    if (!success)
        return false;

    auto allocate_bundle = [&storages, &bundles](const index_t &bi) {
        cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t b;
        const int divx = cslibs_math::common::div<int>(bi[0], 2);
        const int divy = cslibs_math::common::div<int>(bi[1], 2);
        const int modx = cslibs_math::common::mod<int>(bi[0], 2);
        const int mody = cslibs_math::common::mod<int>(bi[1], 2);

        const index_t storage_0_index = {{divx,        divy}};
        const index_t storage_1_index = {{divx + modx, divy}};        /// shifted to the left
        const index_t storage_2_index = {{divx,        divy + mody}}; /// shifted to the bottom
        const index_t storage_3_index = {{divx + modx, divy + mody}}; /// shifted diagonally

        b[0] = storages[0]->get(storage_0_index);
        b[1] = storages[1]->get(storage_1_index);
        b[2] = storages[2]->get(storage_2_index);
        b[3] = storages[3]->get(storage_3_index);
        bundles->insert(bi, b);
    };
    for(const index_t &index : indices)
        allocate_bundle(index);

    map.reset(new cslibs_ndt_2d::static_maps::Gridmap(origin,
                                                      resolution,
                                                      size,
                                                      bundles,
                                                      storages));

    return true;
}
}
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
