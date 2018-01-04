#ifndef CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_3d/static_maps/gridmap.hpp>
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

        using vector_t = std::vector<cslibs_math::statistics::Distribution<3, 3>>;        
        using index_t = std::array<int, 3>;

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
            vector_t storage(i == 0 ? (size[0] * size[1] * size[2]) : ((size[0] + 1) * (size[1] + 1) * (size[2] + 1)));

            for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
                for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                    for (int idz = 0 ; idz < static_cast<int>(rhs->getBundleSize()[2]) ; ++ idz) {
                        index_t bi({idx, idy, idz});
                        if (const typename cslibs_ndt_3d::static_maps::Gridmap::distribution_bundle_t* b =
                                rhs->getDistributionBundle(bi)) {

                            const index_t storage_index = get_storage_index(bi, i);
                            const std::size_t index = (i == 0 ?
                                (storage_index[2] * size[1] + storage_index[1]) * size[0] + storage_index[0] :
                                (storage_index[2] * (size[1] + 1) + storage_index[1]) * (size[0] + 1) + storage_index[0]);

                            storage[index] = *(b->at(i));
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

        const std::array<std::size_t, 3> & size = rhs->getSize();

        using vector_t = std::vector<cslibs_math::statistics::Distribution<3, 3>>;
        using index_t = std::array<int, 3>;

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
            vector_t storage = n[3 + i].as<vector_t>();

            for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
                for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                    for (int idz = 0 ; idz < static_cast<int>(rhs->getBundleSize()[2]) ; ++ idz) {
                        index_t bi({idx, idy, idz});
                        if (const typename cslibs_ndt_3d::static_maps::Gridmap::distribution_bundle_t* b =
                                rhs->getDistributionBundle(bi)) {

                            const index_t storage_index = get_storage_index(bi, i);
                            const std::size_t index = (i == 0) ?
                                ((storage_index[2] * size[1] + storage_index[1]) * size[0] + storage_index[0]) :
                                ((storage_index[2] * (size[1] + 1) + storage_index[1]) * (size[0] + 1) + storage_index[0]);

                            b->at(i)->data() = storage[index];
                        }
                    }
                }
            }
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_3D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
