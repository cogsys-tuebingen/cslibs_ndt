#ifndef CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
#define CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP

#include <cslibs_ndt/common/indexed.hpp>
#include <cslibs_ndt/common/distribution.hpp>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_math/serialization/distribution.hpp>

#include <yaml-cpp/yaml.h>

template<typename std::size_t Size, std::size_t Dim>
YAML::Emitter& operator << (YAML::Emitter& out, const cslibs_ndt::Indexed<cslibs_ndt::Distribution, Size, Dim> &rhs) {
    YAML::Node n;
    n.push_back(rhs.index_);
    n.push_back(rhs.data_.getHandle()->data());
    out << n;
    return out;
}


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
        rhs.data_.data() = n[1].as<typename cslibs_ndt::Distribution<Size>::distribution_t>();

        return true;
    }
};
}

#endif // CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
