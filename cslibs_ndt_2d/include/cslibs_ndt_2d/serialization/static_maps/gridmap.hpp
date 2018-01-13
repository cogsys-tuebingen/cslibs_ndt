#ifndef CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt/common/serialization/storage.hpp>

#include <cslibs_ndt_2d/static_maps/gridmap.hpp>
#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>

#include <yaml-cpp/yaml.h>

namespace YAML {
template <>
struct convert<cslibs_ndt_2d::static_maps::Gridmap::Ptr>
{
    using map_t = cslibs_ndt_2d::static_maps::Gridmap;
    static Node encode(const map_t::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getOrigin());
        n.push_back(rhs->getResolution());

        const std::array<std::size_t, 2> & size = rhs->getSize();
        n.push_back(size);

        using index_t = std::array<int, 2>;
        using distribution_storage_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_ptr_t;

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

            for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
                for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                    const index_t bi({idx, idy});

                    if (const typename map_t::distribution_bundle_t* b = rhs->getDistributionBundle(bi)) {
                        const index_t si = get_storage_index(bi, i);
                        if (!storage->get(si) && b->at(i)->data().getN() > 0)
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
        if (!n.IsSequence() || n.size() != 7)
            return false;

        rhs.reset(new map_t(n[0].as<cslibs_math_2d::Transform2d>(), n[1].as<double>(), n[2].as<std::array<std::size_t, 2>>()));

        using index_t = std::array<int, 2>;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_ptr_t;

        const std::array<std::size_t, 2> & bundle_size = rhs->getBundleSize();
        auto get_bundle_index = [&bundle_size] (const index_t & si) {
            return index_t({{std::max(0, std::min(2 * si[0], static_cast<int>(bundle_size[0] - 1))),
                             std::max(0, std::min(2 * si[1], static_cast<int>(bundle_size[1] - 1)))}});
        };

        for (std::size_t i = 0 ; i < 4 ; ++ i) {
            const distribution_storage_ptr_t & storage = n[3 + i].as<distribution_storage_ptr_t>();

            storage->traverse([&rhs, &i, &get_bundle_index] (const index_t & si, const typename map_t::distribution_t & d) {
                const index_t & bi = get_bundle_index(si);
                if (const typename map_t::distribution_bundle_t* b = rhs->getDistributionBundle(bi))
                    b->at(i)->data() = d.data();
            });
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
