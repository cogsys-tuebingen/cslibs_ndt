#ifndef CSLIBS_NDT_STORAGE_HPP
#define CSLIBS_NDT_STORAGE_HPP

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_ndt/common/indexed.hpp>
#include <yaml-cpp/yaml.h>

namespace cis = cslibs_indexed_storage;

namespace YAML {
template <template <std::size_t> class T, std::size_t Size, std::size_t Dim>
struct convert<std::shared_ptr<cis::Storage<T<Size>,
        std::array<int, Dim>, cis::backend::kdtree::KDTree>>>
{
    static Node encode(const std::shared_ptr<cis::Storage<T<Size>,
                       std::array<int, Dim>, cis::backend::kdtree::KDTree>> &rhs)
    {
        Node n;
        if (!rhs)
            return n;

        using index_t = std::array<int, Dim>;
        using data_t  = T<Size>;
        rhs->traverse([&n](const index_t& index, const data_t& data) {
            cslibs_ndt::Indexed<T, Size, Dim> d(index, data);
            n.push_back(d);
        });

        return n;
    }

    static bool decode(const Node& n, std::shared_ptr<cis::Storage<T<Size>,
                       std::array<int, Dim>, cis::backend::kdtree::KDTree>> &rhs)
    {
        if (!n.IsSequence())
            return false;

        rhs.reset(new cis::Storage<T<Size>, std::array<int, Dim>, cis::backend::kdtree::KDTree>());
        for (std::size_t p = 0 ; p < n.size() ; ++ p) {
            cslibs_ndt::Indexed<T, Size, Dim> d = n[p].as<cslibs_ndt::Indexed<T, Size, Dim>>();
            rhs->insert(d.index_, d.data_);
        }

        return true;
    }
};
}

#endif // CSLIBS_NDT_STORAGE_HPP
