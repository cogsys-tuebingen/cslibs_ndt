#ifndef CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_occupancy_distribution.hpp>
#include <cslibs_ndt/common/serialization/storage.hpp>

#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>

#include <yaml-cpp/yaml.h>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_ndt/common/serialization/filesystem.hpp>
#include <fstream>
#include <thread>
#include <atomic>

namespace cslibs_ndt_2d {
namespace dynamic_maps {
inline bool saveBinary(const cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &map,
                       const std::string &path)
{
    using path_t                     = boost::filesystem::path;
    using paths_t                    = std::array<path_t, 4>;
    using index_t                    = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::index_t;
    using distribution_storage_ptr_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_storage_ptr_t;
    using binary_t                   = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 2, 2>;
    using storages_t                 = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::create_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
    paths_t paths = {{path_root / path_t("store_0.bin"),
                      path_root / path_t("store_1.bin"),
                      path_root / path_t("store_2.bin"),
                      path_root / path_t("store_3.bin")}};

    /// step three: we have our filesystem, now we write out the distributions file by file
    /// meta file
    path_t path_file = path_t("map.yaml");
    {
        std::ofstream out = std::ofstream((path_root / path_file).string());
        YAML::Emitter yaml(out);
        YAML::Node n;
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"] = map->getInitialOrigin();
        n["resolution"] = map->getResolution();
        n["min_index"] = map->getMinDistributionIndex();
        n["max_index"] = map->getMaxDistributionIndex();
        n["bundles"] = indices;
        yaml << n;
    }
    /// step four: write out the storages


    const storages_t storages = {{map->getStorages()[0],
                                  map->getStorages()[1],
                                  map->getStorages()[2],
                                  map->getStorages()[3]}};
    std::array<std::thread, 4> threads;
    std::atomic_bool success(true);
    for(std::size_t i = 0 ; i < 4 ; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::save(storages[i], paths[i]);
        });
    for(std::size_t i = 0 ; i < 4 ; ++i)
        threads[i].join();

    return success;
}

inline bool loadBinary(const std::string &path,
                       cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &map)
{
    using path_t                    = boost::filesystem::path;
    using paths_t                   = std::array<path_t, 4>;
    using index_t                   = std::array<int, 2>;
    using binary_t                  = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 2, 2>;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::check_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
    /// step two: check if the sub folders can be created
    paths_t paths = {{path_root / path_t("store_0.bin"),
                      path_root / path_t("store_1.bin"),
                      path_root / path_t("store_2.bin"),
                      path_root / path_t("store_3.bin")}};

    /// step three: we have our filesystem, now we can load distributions file by file
    /// load meta data
    path_t  path_file = path_t("map.yaml");

    using bundle_storage_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_bundle_storage_t;
    using distribution_storage_array_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;

    std::shared_ptr<bundle_storage_t> bundles(new bundle_storage_t);
    distribution_storage_array_t storages;

    YAML::Node n = YAML::LoadFile((path_root / path_file).string());
    const cslibs_math_2d::Transform2d origin = n["origin"].as<cslibs_math_2d::Transform2d>();
    const double               resolution = n["resolution"].as<double>();
    const index_t              min_index  = n["min_index"].as<index_t>();
    const index_t              max_index  = n["max_index"].as<index_t>();
    const std::vector<index_t> indices    = n["bundles"].as<std::vector<index_t>>();

    std::array<std::thread, 4> threads;
    std::atomic_bool success(true);
    for(std::size_t i = 0 ; i < 4 ; ++i)
        threads[i] = std::thread([&storages, &paths, i, &success](){
            success = success && binary_t::load(paths[i], storages[i]);
        });
    for(std::size_t i = 0 ; i < 4 ; ++i)
        threads[i].join();

    if(!success)
        return false;

    auto allocate_bundle = [&storages, &bundles](const index_t &bi) {
        cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_bundle_t b;
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


    for(const index_t &index : indices) {
        allocate_bundle(index);
    }
    map.reset(new cslibs_ndt_2d::dynamic_maps::OccupancyGridmap(origin, resolution,
                                                       bundles,
                                                       storages,
                                                       min_index,
                                                       max_index));

    return true;
}



inline bool save(const cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &map,
                 const std::string &path)
{
    using path_t  = boost::filesystem::path;
    using paths_t = std::array<path_t, 4>;
    using index_t = std::array<int, 2>;
    using distribution_storage_ptr_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::distribution_storage_ptr_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::create_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
    paths_t paths = {{path_root / path_t("0"), path_root / path_t("1"), path_root / path_t("2"), path_root / path_t("3")}};

    /// step three: we have our filesystem, now we write out the distributions file by file
    /// meta file
    path_t path_file = path_t("map.yaml");
    {
        std::ofstream out = std::ofstream((path_root / path_file).string());
        YAML::Emitter yaml(out);
        YAML::Node n;
        std::vector<index_t> indices;
        map->getBundleIndices(indices);
        n["origin"] = map->getInitialOrigin();
        n["resolution"] = map->getResolution();
        n["min_index"] = map->getMinDistributionIndex();
        n["max_index"] = map->getMaxDistributionIndex();
        n["bundles"] = indices;
        yaml << n;

    }
    for(std::size_t i = 0 ; i < 4 ; ++i) {
        if(!cslibs_ndt::common::serialization::create_directory(paths[i])) {
            return false;
        }
    }

    const distribution_storage_ptr_t storage_0 = map->getStorages()[0];
    const distribution_storage_ptr_t storage_1 = map->getStorages()[1];
    const distribution_storage_ptr_t storage_2 = map->getStorages()[2];
    const distribution_storage_ptr_t storage_3 = map->getStorages()[3];

    cslibs_ndt::save(storage_0, paths[0]);
    cslibs_ndt::save(storage_1, paths[1]);
    cslibs_ndt::save(storage_2, paths[2]);
    cslibs_ndt::save(storage_3, paths[3]);

    return true;
}

inline bool load(cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr &map,
                 const std::string &path)
{
    using path_t  = boost::filesystem::path;
    using paths_t = std::array<path_t, 4>;
    using index_t = std::array<int, 2>;
    using map_t   = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap;
    using distribution_storage_ptr_t   = typename map_t::distribution_storage_ptr_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::check_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
    paths_t paths = {{path_root / path_t("0"), path_root /  path_t("1"), path_root /  path_t("2"), path_root / path_t("3")}};

    /// step three: we have our filesystem, now we can load distributions file by file
    for(std::size_t i = 0 ; i < 4 ; ++i) {
        if(!cslibs_ndt::common::serialization::check_directory(paths[i])) {
            return false;
        }
    }

    /// load meta data
    path_t path_file = path_t("map.yaml");
    index_t min_index, max_index;
    {
        YAML::Node n = YAML::LoadFile((path_root / path_file).string());
        const cslibs_math_2d::Transform2d origin = n["origin"].as<cslibs_math_2d::Transform2d>();
        const double resolution = n["resolution"].as<double>();
        min_index = n["min_index"].as<index_t>();
        max_index = n["max_index"].as<index_t>();
        const std::vector<index_t> indices = n["bundles"].as<std::vector<index_t>>();

        map.reset(new map_t(origin, resolution));

        // allocation stuff, just 4 fun
        for (const index_t& bi : indices)
            map->getDistributionBundle(bi);
    }

    auto get_bundle_index = [&min_index, &max_index] (const index_t & si) {
        return index_t({{std::max(min_index[0], std::min(2 * si[0], max_index[0])),
                         std::max(min_index[1], std::min(2 * si[1], max_index[1]))}});
    };

    for (std::size_t i = 0 ; i < 4 ; ++ i) {
        distribution_storage_ptr_t storage;
        if (!cslibs_ndt::load(storage, paths[i]))
            return false;

        storage->traverse([&map, &i, &get_bundle_index] (const index_t & si, const typename map_t::distribution_t & d) {
            const index_t & bi = get_bundle_index(si);
            if (const typename map_t::distribution_bundle_t* b = map->getDistributionBundle(bi))
                if (b->at(i))
                    *(b->at(i)) = d;
        });
    }

    return true;
}
}
}

namespace YAML {
template <>
struct convert<cslibs_ndt_2d::dynamic_maps::OccupancyGridmap::Ptr>
{
    using map_t = cslibs_ndt_2d::dynamic_maps::OccupancyGridmap;
    static Node encode(const typename map_t::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getInitialOrigin());
        n.push_back(rhs->getResolution());

        using index_t = std::array<int, 2>;
        const index_t min_distribution_index = rhs->getMinDistributionIndex();
        const index_t max_distribution_index = rhs->getMaxDistributionIndex();
        n.push_back(min_distribution_index);
        n.push_back(max_distribution_index);

        using distribution_storage_t = typename map_t::distribution_storage_t;
        using distribution_storage_ptr_t = typename map_t::distribution_storage_ptr_t;

        auto divx = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[0], 2); };
        auto divy = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[1], 2); };
        auto modx = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[0], 2); };
        auto mody = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[1], 2); };

        auto get_storage_index = [&divx, &divy, &modx, &mody](const index_t & bi, const std::size_t i) {
            return index_t({{(i % 2 == 0) ? divx(bi) : (divx(bi) + modx(bi)),
                             (i < 2) ? divy(bi) : (divy(bi) + mody(bi))}});
        };

        for (std::size_t i = 0 ; i < 4 ; ++ i) {
            const distribution_storage_ptr_t storage(new distribution_storage_t());

            for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
                for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
                    const index_t bi({idx, idy});

                    if (const typename map_t::distribution_bundle_t* b = rhs->getDistributionBundle(bi)) {
                        const index_t si = get_storage_index(bi, i);
                        if (!storage->get(si) && b->at(i) && (b->at(i)->numFree() > 0 || b->at(i)->numOccupied() > 0))
                            storage->insert(si, *(b->at(i)));
                    }
                }
            }

            n.push_back(storage);
        }

        return n;
    }

    static bool decode(const Node& n, typename map_t::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 8)
            return false;

        rhs.reset(new map_t(n[0].as<cslibs_math_2d::Transform2d>(), n[1].as<double>()));

        using index_t = std::array<int, 2>;
        const index_t min_distribution_index = n[2].as<index_t>();
        const index_t max_distribution_index = n[3].as<index_t>();

        using distribution_storage_ptr_t = typename map_t::distribution_storage_ptr_t;

        auto get_bundle_index = [&min_distribution_index, &max_distribution_index] (const index_t & si) {
            return index_t({{std::max(min_distribution_index[0], std::min(2 * si[0], max_distribution_index[0])),
                             std::max(min_distribution_index[1], std::min(2 * si[1], max_distribution_index[1]))}});
        };

        for (std::size_t i = 0 ; i < 4 ; ++ i) {
            const distribution_storage_ptr_t & storage = n[4 + i].as<distribution_storage_ptr_t>();

            storage->traverse([&rhs, &i, &get_bundle_index] (const index_t & si, const typename map_t::distribution_t & d) {
                const index_t & bi = get_bundle_index(si);
                if (typename map_t::distribution_bundle_t* b = rhs->getDistributionBundle(bi))
                    if (b->at(i))
                        *(b->at(i)) = d;
            });
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
