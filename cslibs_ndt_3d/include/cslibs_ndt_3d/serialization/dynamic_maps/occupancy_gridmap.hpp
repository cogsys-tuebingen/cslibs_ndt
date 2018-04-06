#ifndef CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_ndt/serialization/filesystem.hpp>
#include <cslibs_ndt/serialization/storage.hpp>
#include <cslibs_ndt/serialization/indexed_occupancy_distribution.hpp>

#include <cslibs_math_3d/serialization/transform.hpp>
#include <cslibs_math/serialization/array.hpp>

#include <yaml-cpp/yaml.h>

#include <fstream>
#include <thread>
#include <atomic>

namespace cslibs_ndt_3d {
namespace dynamic_maps {
inline bool saveBinary(const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &map,
                       const std::string &path)
{
    using path_t     = boost::filesystem::path;
    using paths_t    = std::array<path_t, 8>;
    using index_t    = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::index_t;
    using storages_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;
    using binary_t   = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 3, 3>;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::create_directory(path_root))
        return false;

    /// step two: identity subfolders
    const paths_t paths = {{path_root / path_t("store_0.bin"),
                            path_root / path_t("store_1.bin"),
                            path_root / path_t("store_2.bin"),
                            path_root / path_t("store_3.bin"),
                            path_root / path_t("store_4.bin"),
                            path_root / path_t("store_5.bin"),
                            path_root / path_t("store_6.bin"),
                            path_root / path_t("store_7.bin")}};

    /// step three: we have our filesystem, now we write out the distributions file by file
    /// meta file
    const path_t path_file = path_t("map.yaml");
    {
        std::ofstream out((path_root / path_file).string(), std::fstream::trunc);
        YAML::Emitter yaml(out);
        YAML::Node n;
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"]     = map->getInitialOrigin();
        n["resolution"] = map->getResolution();
        n["min_index"]  = map->getMinDistributionIndex();
        n["max_index"]  = map->getMaxDistributionIndex();
        n["bundles"]    = indices;
        yaml << n;
    }

    /// step four: write out the storages
    const storages_t storages = {{map->getStorages()[0],
                                  map->getStorages()[1],
                                  map->getStorages()[2],
                                  map->getStorages()[3],
                                  map->getStorages()[4],
                                  map->getStorages()[5],
                                  map->getStorages()[6],
                                  map->getStorages()[7]}};

    std::array<std::thread, 8> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < 8 ; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::save(storages[i], paths[i]);
        });
    for (std::size_t i = 0 ; i < 8 ; ++i)
        threads[i].join();

    return success;
}

inline bool loadBinary(const std::string &path,
                       cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &map)
{
    using path_t           = boost::filesystem::path;
    using paths_t          = std::array<path_t, 8>;
    using index_t          = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::index_t;
    using binary_t         = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 3, 3>;
    using bundle_storage_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_bundle_storage_t;
    using storages_t       = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if (!cslibs_ndt::common::serialization::check_directory(path_root))
        return false;

    /// step two: identity subfolders
    const paths_t paths = {{path_root / path_t("store_0.bin"),
                            path_root / path_t("store_1.bin"),
                            path_root / path_t("store_2.bin"),
                            path_root / path_t("store_3.bin"),
                            path_root / path_t("store_4.bin"),
                            path_root / path_t("store_5.bin"),
                            path_root / path_t("store_6.bin"),
                            path_root / path_t("store_7.bin")}};

    /// step three: we have our filesystem, now we can load distributions file by file
    for (std::size_t i = 0 ; i < 8 ; ++i)
        if (!cslibs_ndt::common::serialization::check_file(paths[i]))
            return false;

    /// load meta data
    path_t  path_file = path_t("map.yaml");

    std::shared_ptr<bundle_storage_t> bundles(new bundle_storage_t);
    storages_t storages;

    YAML::Node n = YAML::LoadFile((path_root / path_file).string());
    const cslibs_math_3d::Transform3d origin     = n["origin"].as<cslibs_math_3d::Transform3d>();
    const double                      resolution = n["resolution"].as<double>();
    const index_t                     min_index  = n["min_index"].as<index_t>();
    const index_t                     max_index  = n["max_index"].as<index_t>();
    const std::vector<index_t>        indices    = n["bundles"].as<std::vector<index_t>>();

    std::array<std::thread, 8> threads;
    std::atomic_bool success(true);
    for (std::size_t i = 0 ; i < 8 ; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::load(paths[i], storages[i]);
        });
    for (std::size_t i = 0 ; i < 8 ; ++i)
        threads[i].join();

    if (!success)
        return false;

    auto allocate_bundle = [&storages, &bundles](const index_t &bi) {
        cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_bundle_t b;
        const int divx = cslibs_math::common::div<int>(bi[0], 2);
        const int divy = cslibs_math::common::div<int>(bi[1], 2);
        const int divz = cslibs_math::common::div<int>(bi[2], 2);
        const int modx = cslibs_math::common::mod<int>(bi[0], 2);
        const int mody = cslibs_math::common::mod<int>(bi[1], 2);
        const int modz = cslibs_math::common::mod<int>(bi[2], 2);

        const index_t storage_0_index = {{divx,        divy,        divz}};
        const index_t storage_1_index = {{divx + modx, divy,        divz}};
        const index_t storage_2_index = {{divx,        divy + mody, divz}};
        const index_t storage_3_index = {{divx + modx, divy + mody, divz}};
        const index_t storage_4_index = {{divx,        divy,        divz + modz}};
        const index_t storage_5_index = {{divx + modx, divy,        divz + modz}};
        const index_t storage_6_index = {{divx,        divy + mody, divz + modz}};
        const index_t storage_7_index = {{divx + modx, divy + mody, divz + modz}};

        b[0] = storages[0]->get(storage_0_index);
        b[1] = storages[1]->get(storage_1_index);
        b[2] = storages[2]->get(storage_2_index);
        b[3] = storages[3]->get(storage_3_index);
        b[4] = storages[4]->get(storage_4_index);
        b[5] = storages[5]->get(storage_5_index);
        b[6] = storages[6]->get(storage_6_index);
        b[7] = storages[7]->get(storage_7_index);
        bundles->insert(bi, b);
    };
    for (const index_t &index : indices)
        allocate_bundle(index);

    map.reset(new cslibs_ndt_3d::dynamic_maps::OccupancyGridmap(origin,
                                                                resolution,
                                                                min_index,
                                                                max_index,
                                                                bundles,
                                                                storages));

    return true;
}
}
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
