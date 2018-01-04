#ifndef CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>
#include <yaml-cpp/yaml.h>

namespace YAML {
template <>
struct convert<cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr>
{
    static Node encode(const cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getInitialOrigin());
        n.push_back(rhs->getResolution());

        using index_t = std::array<int, 2>;
        index_t min_distribution_index = rhs->getMinDistributionIndex();
        index_t max_distribution_index = rhs->getMaxDistributionIndex();
        n.push_back(min_distribution_index);
        n.push_back(max_distribution_index);

        using distribution_storage_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_ptr_t;
        using distribution_storage_array_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_array_t;
        distribution_storage_array_t storage({{distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t)}});

        for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
            for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
                index_t bi({idx, idy});
                if (const typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi)) {

                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    for (std::size_t i = 0 ; i < 4 ; ++ i)
                        if (!storage[i]->get(storage_indices[i]) && b->at(i))
                            storage[i]->insert(storage_indices[i], *(b->at(i)));
                }
            }
        }

        for (std::size_t i = 0 ; i < 4 ; ++ i)
            n.push_back(storage[i]);

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt_2d::dynamic_maps::Gridmap::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 8)
            return false;

        rhs.reset(new cslibs_ndt_2d::dynamic_maps::Gridmap(
                      n[0].as<cslibs_math_2d::Transform2d>(), n[1].as<double>()));

        using index_t = std::array<int, 2>;
        index_t min_distribution_index = n[2].as<index_t>();
        index_t max_distribution_index = n[3].as<index_t>();

        using distribution_storage_array_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_array_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_storage_ptr_t;
        distribution_storage_array_t storage({n[4].as<distribution_storage_ptr_t>(),
                                              n[5].as<distribution_storage_ptr_t>(),
                                              n[6].as<distribution_storage_ptr_t>(),
                                              n[7].as<distribution_storage_ptr_t>()});

        for (int idx = min_distribution_index[0] ; idx <= max_distribution_index[0] ; ++ idx) {
            for (int idy = min_distribution_index[1] ; idy <= max_distribution_index[1] ; ++ idy) {
                index_t bi({idx, idy});
                if (typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi)) {

                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    for (std::size_t i = 0 ; i < 4 ; ++ i)
                        if (storage[i]->get(storage_indices[i]) && b->at(i))
                            b->at(i)->data() = storage[i]->get(storage_indices[i])->data();
                }
            }
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_DYNAMIC_MAPS_GRIDMAP_HPP
