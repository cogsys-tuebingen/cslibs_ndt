#ifndef CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_occupancy_distribution.hpp>
#include <cslibs_ndt/common/serialization/storage.hpp>

#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_math_3d/serialization/transform.hpp>

#include <yaml-cpp/yaml.h>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_ndt/common/serialization/filesystem.hpp>
#include <fstream>

namespace cslibs_ndt_3d {
namespace dynamic_maps {
inline bool saveBinary(const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &map,
                       const std::string &path)
{
    using path_t                        = boost::filesystem::path;
    using paths_t                       = std::array<path_t, 8>;
    using index_t                       = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::index_t;
    using distribution_storage_array_t  = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;
    using binary_t                      = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 3, 3>;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::create_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
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

    const distribution_storage_array_t storages = {{map->getStorages()[0],
                                                    map->getStorages()[1],
                                                    map->getStorages()[2],
                                                    map->getStorages()[3],
                                                    map->getStorages()[4],
                                                    map->getStorages()[5],
                                                    map->getStorages()[6],
                                                    map->getStorages()[7]}};

    for(std::size_t i = 0 ; i < 8 ; ++i) {
        if(!binary_t::save(storages[i], paths[i]))
            return false;
    }
    return true;
}

inline bool loadBinary(const std::string &path,
                       cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &map)
{
    using path_t                        = boost::filesystem::path;
    using paths_t                       = std::array<path_t, 8>;
    using index_t                       = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::index_t;
    using binary_t                      = cslibs_ndt::binary<cslibs_ndt::OccupancyDistribution, 3, 3>;
    using bundle_storage_t              = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_bundle_storage_t;
    using distribution_storage_array_t  = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_storage_array_t;

    /// step one: check if the root diretory exists
    path_t path_root(path);
    if(!cslibs_ndt::common::serialization::check_directory(path_root)) {
        return false;
    }

    /// step two: check if the sub folders can be created
    /// step two: check if the sub folders can be created
    const paths_t paths = {{path_root / path_t("store_0.bin"),
                            path_root / path_t("store_1.bin"),
                            path_root / path_t("store_2.bin"),
                            path_root / path_t("store_3.bin"),
                            path_root / path_t("store_4.bin"),
                            path_root / path_t("store_5.bin"),
                            path_root / path_t("store_6.bin"),
                            path_root / path_t("store_7.bin")}};

    /// step three: we have our filesystem, now we can load distributions file by file
    for(std::size_t i = 0 ; i < 4 ; ++i) {
        if(!cslibs_ndt::common::serialization::check_file(paths[i])) {
            return false;
        }
    }

    /// load meta data
    path_t  path_file = path_t("map.yaml");

    std::shared_ptr<bundle_storage_t> bundles(new bundle_storage_t);
    distribution_storage_array_t storages;

    YAML::Node n = YAML::LoadFile((path_root / path_file).string());
    const cslibs_math_3d::Transform3d origin = n["origin"].as<cslibs_math_3d::Transform3d>();
    const double                      resolution = n["resolution"].as<double>();
    const index_t                     min_index  = n["min_index"].as<index_t>();
    const index_t                     max_index  = n["max_index"].as<index_t>();
    const std::vector<index_t>        indices    = n["bundles"].as<std::vector<index_t>>();

    for(int i = 0 ; i < 8 ; ++i) {
        if (!binary_t::load(paths[i], storages[i]))
            return false;
    }

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


    for(const index_t &index : indices) {
        allocate_bundle(index);
    }
    map.reset(new cslibs_ndt_3d::dynamic_maps::OccupancyGridmap(origin, resolution,
                                                                bundles,
                                                                storages,
                                                                min_index,
                                                                max_index));

    return true;
}
}
}

namespace YAML {
template <>
struct convert<cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr>
{
    using map_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap;
    static Node encode(const typename map_t::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getInitialOrigin());
        n.push_back(rhs->getResolution());

        using index_t = std::array<int, 3>;
        const index_t min_distribution_index = rhs->getMinDistributionIndex();
        const index_t max_distribution_index = rhs->getMaxDistributionIndex();
        n.push_back(min_distribution_index);
        n.push_back(max_distribution_index);

        using distribution_storage_t = typename map_t::distribution_storage_t;
        using distribution_storage_ptr_t = typename map_t::distribution_storage_ptr_t;

        auto divx = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[0], 2); };
        auto divy = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[1], 2); };
        auto divz = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[2], 2); };
        auto modx = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[0], 2); };
        auto mody = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[1], 2); };
        auto modz = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[2], 2); };

        auto get_storage_index = [&divx, &divy, &divz, &modx, &mody, &modz] (const index_t & bi, const std::size_t i) {
            return index_t({{(i % 2 == 0) ? divx(bi) : (divx(bi) + modx(bi)),
                             ((i / 2) % 2 == 0) ? divy(bi) : (divy(bi) + mody(bi)),
                             (i < 4) ? divz(bi) : (divz(bi) + modz(bi))}});
        };

        for (std::size_t i = 0 ; i < 8 ; ++ i) {
            const distribution_storage_ptr_t storage(new distribution_storage_t());

            for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
                for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
                    for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz) {
                        const index_t bi({idx, idy, idz});

                        if (const typename map_t::distribution_bundle_t* b = rhs->getDistributionBundle(bi)) {
                            const index_t si = get_storage_index(bi, i);
                            if (!storage->get(si) && b->at(i) && (b->at(i)->numFree() > 0 || b->at(i)->numOccupied() > 0))
                                storage->insert(si, *(b->at(i)));
                        }
                    }
                }
            }

            n.push_back(storage);
        }

        return n;
    }

    static bool decode(const Node& n, typename map_t::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 12)
            return false;

        rhs.reset(new map_t(n[0].as<cslibs_math_3d::Transform3d>(), n[1].as<double>()));

        using index_t = std::array<int, 3>;
        const index_t min_distribution_index = n[2].as<index_t>();
        const index_t max_distribution_index = n[3].as<index_t>();

        using distribution_storage_ptr_t = typename map_t::distribution_storage_ptr_t;

        auto get_bundle_index = [&min_distribution_index, &max_distribution_index] (const index_t & si) {
            return index_t({{std::max(min_distribution_index[0], std::min(2 * si[0], max_distribution_index[0])),
                             std::max(min_distribution_index[1], std::min(2 * si[1], max_distribution_index[1])),
                             std::max(min_distribution_index[2], std::min(2 * si[2], max_distribution_index[2]))}});
        };

        for (std::size_t i = 0 ; i < 8 ; ++ i) {
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

#endif // CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
