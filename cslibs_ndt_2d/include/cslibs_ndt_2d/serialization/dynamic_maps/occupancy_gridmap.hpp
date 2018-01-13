#ifndef CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_occupancy_distribution.hpp>
#include <cslibs_ndt/common/serialization/storage.hpp>

#include <cslibs_ndt_2d/dynamic_maps/occupancy_gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>

#include <yaml-cpp/yaml.h>

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
