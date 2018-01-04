#ifndef CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_2d/static_maps/gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>
#include <yaml-cpp/yaml.h>

namespace YAML {
template <>
struct convert<cslibs_ndt_2d::static_maps::Gridmap::Ptr>
{
    static Node encode(const cslibs_ndt_2d::static_maps::Gridmap::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getOrigin());
        n.push_back(rhs->getResolution());
        n.push_back(rhs->getSize());

        using distribution_storage_t =
        typename cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_ptr_t;
        using distribution_storage_array_t =
        typename cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_array_t;
        distribution_storage_array_t storage({{distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t),
                                               distribution_storage_ptr_t(new distribution_storage_t)}});

        using index_t = std::array<int, 2>;
        for (int idx = 0 ; idx <= rhs->getBundleSize()[0] ; ++ idx) {
            for (int idy = 0 ; idy <= rhs->getBundleSize()[1] ; ++ idy) {
                index_t bi({idx, idy});
                if (const typename cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t* d =
                        rhs->getDistributionBundle(bi)) {
                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    for (std::size_t i = 0 ; i < 4 ; ++ i)
                        if (!storage[i]->get(storage_indices[i]) && d->at(i))
                            storage[i]->insert(storage_indices[i], *(d->at(i)));
                }
            }
        }

        for (std::size_t i = 0 ; i < 4 ; ++ i)
            n.push_back(storage[i]);

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt_2d::static_maps::Gridmap::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 7)
            return false;

        rhs.reset(new cslibs_ndt_2d::static_maps::Gridmap(
                      n[0].as<cslibs_math_2d::Transform2d>(), n[1].as<double>(), n[2].as<std::size_t>()));

        using distribution_storage_array_t =
        typename cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_array_t;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_2d::static_maps::Gridmap::distribution_storage_ptr_t;
        distribution_storage_array_t storage({n[3].as<distribution_storage_ptr_t>(),
                                              n[4].as<distribution_storage_ptr_t>(),
                                              n[5].as<distribution_storage_ptr_t>(),
                                              n[6].as<distribution_storage_ptr_t>()});

        using index_t = std::array<int, 2>;
        for (int idx = 0 ; idx <= rhs->getBundleSize()[0] ; ++ idx) {
            for (int idy = 0 ; idy <= rhs->getBundleSize()[1] ; ++ idy) {
                index_t bi({idx, idy});
                if (typename cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi)) {
                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    for (std::size_t i = 0 ; i < 4 ; ++ i)
                        b->at(i) = storage[i]->get(storage_indices[i]);
                }
            }
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
