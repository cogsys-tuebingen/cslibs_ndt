#ifndef CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
#define CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP

#include <cslibs_ndt/common/distribution.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_math/serialization/distribution.hpp>
#include <yaml-cpp/yaml.h>

namespace cslibs_ndt {
template <typename std::size_t Size, std::size_t Dim>
struct IndexedDistribution {
    std::array<int, Dim> index_;
    Distribution<Size>   distribution_;

    IndexedDistribution() = default;

    IndexedDistribution(
            const std::array<int, Dim> & index,
            const Distribution<Size>   & distribution) :
        index_(index),
        distribution_(distribution)
    {
    }
};
}

namespace YAML {
template <typename std::size_t Size, std::size_t Dim>
struct convert<cslibs_ndt::IndexedDistribution<Size, Dim>>
{
    static Node encode(const cslibs_ndt::IndexedDistribution<Size, Dim> &rhs)
    {
        Node n;

        n.push_back(rhs.index_);
        n.push_back(rhs.distribution_.getHandle()->data());

        return n;
    }

    static bool decode(const Node& n, cslibs_ndt::IndexedDistribution<Size, Dim> &rhs)
    {
        if (!n.IsSequence() || n.size() != 2)
            return false;

        rhs.index_               = n[0].as<std::array<int, Dim>>();
        rhs.distribution_.data() = n[1].as<cslibs_math::statistics::Distribution<Size, 3>>();

        return true;
    }
};
}

namespace cis = cslibs_indexed_storage;

namespace YAML {
template <std::size_t Size, std::size_t Dim, template<typename, typename, typename...> class backend_t>
struct convert<std::shared_ptr<cis::Storage<cslibs_ndt::Distribution<Size>, std::array<int, Dim>, backend_t>>>
{
    static Node encode(const std::shared_ptr<cis::Storage<cslibs_ndt::Distribution<Size>, std::array<int, Dim>, backend_t>> &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        using index_t = std::array<int, Dim>;
        using data_t  = cslibs_ndt::Distribution<Size>;
        rhs->traverse([&n](const index_t& index, const data_t& data) {
            cslibs_ndt::IndexedDistribution<Size, Dim> d(index, data);
            n.push_back(d);
        });

        return n;
    }

    static bool decode(const Node& n, std::shared_ptr<cis::Storage<cslibs_ndt::Distribution<Size>, std::array<int, Dim>, backend_t>> &rhs)
    {
        if (!n.IsSequence())
            return false;

        rhs.reset(new cis::Storage<cslibs_ndt::Distribution<Size>, std::array<int, Dim>, backend_t>());
        for (std::size_t p = 0 ; p < n.size() ; ++ p) {
            cslibs_ndt::IndexedDistribution<Size, Dim> d = n[p].as<cslibs_ndt::IndexedDistribution<Size, Dim>>();
            rhs->insert(d.index_, d.distribution_);
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_SERIALIZATION_INDEXED_DISTRIBUTION_HPP
