#ifndef CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP

#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>
#include <cslibs_ndt_2d/static_maps/gridmap.hpp>
#include <cslibs_math_2d/serialization/transform.hpp>

#include <eigen3/Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(cslibs_math::statistics::Distribution<2, 3>)

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

        const std::array<std::size_t, 2> & size = rhs->getSize();
        n.push_back(size);

        using vector_t = std::vector<cslibs_math::statistics::Distribution<2, 3>>;
        std::array<vector_t, 4> storage = {vector_t(size[0] * size[1]),
                                           vector_t((size[0] + 1) * (size[1] + 1)),
                                           vector_t((size[0] + 1) * (size[1] + 1)),
                                           vector_t((size[0] + 1) * (size[1] + 1))};

        using index_t = std::array<int, 2>;
        for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
            for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                index_t bi({idx, idy});
                if (const typename cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi)) {

                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    storage[0][storage_indices[0][1] * size[0] + storage_indices[0][0]] = *(b->at(0));
                    for (std::size_t i = 1 ; i < 4 ; ++ i)
                        storage[i][storage_indices[i][1] * (size[0] + 1) + storage_indices[i][0]] = *(b->at(i));
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
                      n[0].as<cslibs_math_2d::Transform2d>(), n[1].as<double>(), n[2].as<std::array<std::size_t, 2>>()));

        const std::array<std::size_t, 2> & size = rhs->getSize();
        using vector_t = std::vector<cslibs_math::statistics::Distribution<2, 3>>;
        std::array<vector_t, 4> storage = {n[3].as<vector_t>(), n[4].as<vector_t>(), n[5].as<vector_t>(), n[6].as<vector_t>()};

        using index_t = std::array<int, 2>;
        for (int idx = 0 ; idx < static_cast<int>(rhs->getBundleSize()[0]) ; ++ idx) {
            for (int idy = 0 ; idy < static_cast<int>(rhs->getBundleSize()[1]) ; ++ idy) {
                index_t bi({idx, idy});
                if (typename cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t* b =
                        rhs->getDistributionBundle(bi)) {

                    const int divx = cslibs_math::common::div<int>(bi[0], 2);
                    const int divy = cslibs_math::common::div<int>(bi[1], 2);
                    const int modx = cslibs_math::common::mod<int>(bi[0], 2);
                    const int mody = cslibs_math::common::mod<int>(bi[1], 2);

                    const std::array<index_t, 4> storage_indices =
                    {{{divx, divy}, {divx + modx, divy}, {divx, divy + mody}, {divx + modx, divy + mody}}};

                    b->at(0)->data() = storage[0][storage_indices[0][1] * size[0] + storage_indices[0][0]];
                    for (std::size_t i = 1 ; i < 4 ; ++ i)
                        b->at(i)->data() = storage[i][storage_indices[i][1] * (size[0] + 1) + storage_indices[i][0]];
               }
            }
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_2D_SERIALIZATION_STATIC_MAPS_GRIDMAP_HPP
