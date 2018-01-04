#ifndef CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_math_3d/serialization/transform.hpp>
#include <yaml-cpp/yaml.h>

namespace YAML {
template <>
struct convert<cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr>
{
    static Node encode(const cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr &rhs)
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

        using distribution_storage_t =
        typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_storage_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_storage_ptr_t;

        auto divx = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[0], 2); };
        auto divy = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[1], 2); };
        auto divz = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[2], 2); };
        auto modx = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[0], 2); };
        auto mody = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[1], 2); };
        auto modz = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[2], 2); };

        auto get_storage_index = [&divx, &divy, &divz, &modx, &mody, &modz](const index_t & bi, const std::size_t i) {
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

                        if (const typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_bundle_t* b =
                                rhs->getDistributionBundle(bi))
                            if (b->at(i)->data().getN() > 0)
                                storage->insert(get_storage_index(bi, i), *(b->at(i)));
                    }
                }
            }

            n.push_back(storage);
        }

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 12)
            return false;

        rhs.reset(new cslibs_ndt_3d::dynamic_maps::Gridmap(
                      n[0].as<cslibs_math_3d::Transform3d>(), n[1].as<double>()));

        using index_t = std::array<int, 3>;
        const index_t min_distribution_index = n[2].as<index_t>();
        const index_t max_distribution_index = n[3].as<index_t>();

        using distribution_storage_ptr_t =
        typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_storage_ptr_t;

        auto divx = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[0], 2); };
        auto divy = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[1], 2); };
        auto divz = [](const index_t & bi) { return cslibs_math::common::div<int>(bi[2], 2); };
        auto modx = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[0], 2); };
        auto mody = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[1], 2); };
        auto modz = [](const index_t & bi) { return cslibs_math::common::mod<int>(bi[2], 2); };

        auto get_storage_index = [&divx, &divy, &divz, &modx, &mody, &modz](const index_t & bi, const std::size_t i) {
            return index_t({{(i % 2 == 0) ? divx(bi) : (divx(bi) + modx(bi)),
                             ((i / 2) % 2 == 0) ? divy(bi) : (divy(bi) + mody(bi)),
                             (i < 4) ? divz(bi) : (divz(bi) + modz(bi))}});
        };

        for (std::size_t i = 0 ; i < 8 ; ++ i) {
            const distribution_storage_ptr_t & storage = n[4 + i].as<distribution_storage_ptr_t>();

            for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
                for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
                    for (int idz = min_distribution_index[2] ; idz <= max_distribution_index[2] ; ++ idz) {
                        const index_t bi({idx, idy, idz});

                        if (const typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_bundle_t* b =
                                rhs->getDistributionBundle(bi)) {
                            const index_t index = get_storage_index(bi, i);
                            if (auto d = storage->get(index))
                                if (d->data().getN() > 0)
                                    b->at(i)->data() = d->data();
                        }
                    }
                }
            }
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
