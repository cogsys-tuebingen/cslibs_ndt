#ifndef CSLIBS_NDT_SERIALIZATION_INDEXED_OCCUPANCY_DISTRIBUTION_HPP
#define CSLIBS_NDT_SERIALIZATION_INDEXED_OCCUPANCY_DISTRIBUTION_HPP

#include <cslibs_ndt/common/indexed.hpp>
#include <cslibs_ndt/common/occupancy_distribution.hpp>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_math/serialization/distribution.hpp>

#include <yaml-cpp/yaml.h>

namespace YAML {
template <typename std::size_t Size, std::size_t Dim>
struct convert<cslibs_ndt::Indexed<cslibs_ndt::OccupancyDistribution, Size, Dim>>
{
    static Node encode(const cslibs_ndt::Indexed<cslibs_ndt::OccupancyDistribution, Size, Dim> &rhs)
    {
        Node n;

        n.push_back(rhs.index_);
        n.push_back(rhs.data_.numFree());
        n.push_back(rhs.data_.numOccupied());
        if (rhs.data_.getDistribution())
            n.push_back(*(rhs.data_.getDistribution()));

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt::Indexed<cslibs_ndt::OccupancyDistribution, Size, Dim> &rhs)
    {
        if (!n.IsSequence() || n.size() < 3 || n.size() > 4)
            return false;

        rhs.index_ = n[0].as<std::array<int, Dim>>();
        if (n.size() == 3)
            rhs.data_ = cslibs_ndt::OccupancyDistribution<Size>(n[1].as<std::size_t>());
        else
            rhs.data_ = cslibs_ndt::OccupancyDistribution<Size>(n[1].as<std::size_t>(),
                    n[3].as<typename cslibs_ndt::OccupancyDistribution<Size>::distribution_t>());

        return true;
    }
};
}

#endif // CSLIBS_NDT_SERIALIZATION_INDEXED_OCCUPANCY_DISTRIBUTION_HPP
