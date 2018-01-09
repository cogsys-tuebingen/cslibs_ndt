#ifndef CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
#define CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP

#include <cslibs_ndt/common/indexed.hpp>
#include <cslibs_ndt/common/distribution.hpp>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_math/serialization/distribution.hpp>

#include <yaml-cpp/yaml.h>

namespace YAML {
template <typename std::size_t Size, std::size_t Dim>
struct convert<cslibs_ndt::Indexed<cslibs_ndt::Distribution, Size, Dim>>
{
    static Node encode(const cslibs_ndt::Indexed<cslibs_ndt::Distribution, Size, Dim> &rhs)
    {
        Node n;

        n.push_back(rhs.index_);
        n.push_back(rhs.data_.getHandle()->data());

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt::Indexed<cslibs_ndt::Distribution, Size, Dim> &rhs)
    {
        if (!n.IsSequence() || n.size() != 2)
            return false;

        rhs.index_       = n[0].as<std::array<int, Dim>>();
        rhs.data_.data() = n[1].as<cslibs_math::statistics::Distribution<Size, 3>>();

        return true;
    }
};
}

#endif // CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
