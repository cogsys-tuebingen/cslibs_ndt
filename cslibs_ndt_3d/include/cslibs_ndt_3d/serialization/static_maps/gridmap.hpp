#ifndef CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_3d/static_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_math_3d/serialization/transform.hpp>

#include <eigen3/Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(cslibs_math::statistics::Distribution<3, 3>)

#include <yaml-cpp/yaml.h>

namespace YAML {
template <>
struct convert<cslibs_ndt_3d::static_maps::Gridmap::Ptr>
{
    static Node encode(const cslibs_ndt_3d::static_maps::Gridmap::Ptr &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        n.push_back(rhs->getOrigin());
        n.push_back(rhs->getResolution());

        const std::array<std::size_t, 3> & size = rhs->getSize();
        n.push_back(size);

        using index_t = std::array<int, 3>;
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

        auto get_storage_index = [&divx, &divy, &divz, &modx, &mody, &modz] (const index_t & bi, const std::size_t i) {
            return index_t({{(i % 2 == 0) ? divx(bi) : (divx(bi) + modx(bi)),
                             ((i / 2) % 2 == 0) ? divy(bi) : (divy(bi) + mody(bi)),
                             (i < 4) ? divz(bi) : (divz(bi) + modz(bi))}});
        };

        for (std::size_t i = 0 ; i < 8 ; ++ i) {
            const distribution_storage_ptr_t storage(new distribution_storage_t());

            for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
                for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                    for (int idz = 0 ; idz < static_cast<int>(rhs->getBundleSize()[2]) ; ++ idz) {
                        const index_t bi({idx, idy, idz});

                        if (const typename cslibs_ndt_3d::static_maps::Gridmap::distribution_bundle_t* b =
                                rhs->getDistributionBundle(bi)) {
                            const index_t index = get_storage_index(bi, i);
                            if (b->at(i)->data().getN() > 0 && !storage->get(index))
                                storage->insert(index, *(b->at(i)));
                        }
                    }
                }
            }

            n.push_back(storage);
        }

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt_3d::static_maps::Gridmap::Ptr &rhs)
    {
        if (!n.IsSequence() || n.size() != 11)
            return false;

        rhs.reset(new cslibs_ndt_3d::static_maps::Gridmap(
                      n[0].as<cslibs_math_3d::Transform3d>(), n[1].as<double>(), n[2].as<std::array<std::size_t, 3>>()));

        using index_t = std::array<int, 3>;
        using distribution_storage_ptr_t =
        typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_storage_ptr_t;

        const std::array<std::size_t, 3> & bundle_size = rhs->getBundleSize();
        auto get_bundle_index = [&bundle_size] (const index_t & si, const std::size_t & i) {
            return index_t({{std::max(0, std::min(2 * si[0], static_cast<int>(bundle_size[0] - 1))),
                             std::max(0, std::min(2 * si[1], static_cast<int>(bundle_size[1] - 1))),
                             std::max(0, std::min(2 * si[2], static_cast<int>(bundle_size[2] - 1)))}});
        };

        for (std::size_t i = 0 ; i < 8 ; ++ i) {
            const distribution_storage_ptr_t & storage = n[3 + i].as<distribution_storage_ptr_t>();

            storage->traverse([&rhs, &i, &get_bundle_index] (const index_t & si, const cslibs_ndt::Distribution<3> & d) {
                const index_t & bi = get_bundle_index(si, i);
                if (const typename cslibs_ndt_3d::static_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi))
                    b->at(i)->data() = d.data();
            });
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
