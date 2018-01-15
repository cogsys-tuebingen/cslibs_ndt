#ifndef CSLIBS_NDT_STORAGE_HPP
#define CSLIBS_NDT_STORAGE_HPP

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_ndt/common/indexed.hpp>
#include <cslibs_ndt/common/serialization/filesystem.hpp>

#include <fstream>
#include <yaml-cpp/yaml.h>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
template <template <std::size_t> class T, std::size_t Size, std::size_t Dim>
inline void save(const std::shared_ptr<cis::Storage<T<Size>, std::array<int, Dim>,
                 cis::backend::kdtree::KDTree>> &storage,
                 const boost::filesystem::path &directory)
{
    using index_t = std::array<int, Dim>;
    using data_t  = T<Size>;

    auto write = [&directory] (const index_t &index, const data_t &data) {
        std::string path = "d";
        for (std::size_t i = 0 ; i < Dim ; ++ i) {
            path += std::to_string(index[i]);
            if (i < Dim - 1)
                path += "_";
        }
        boost::filesystem::path file =  directory / boost::filesystem::path(path);

        cslibs_ndt::Indexed<T, Size, Dim> indexed(index, data);
        std::ofstream out(file.string());
        out << YAML::Node(indexed);
    };
    storage->traverse(write);
}

template <template <std::size_t> class T, std::size_t Size, std::size_t Dim>
inline bool load(std::shared_ptr<cis::Storage<T<Size>, std::array<int, Dim>,
                 cis::backend::kdtree::KDTree>> &storage,
                 const boost::filesystem::path &directory)
{
    storage.reset(new cis::Storage<T<Size>, std::array<int, Dim>, cis::backend::kdtree::KDTree>());

    for (boost::filesystem::directory_entry & entry : boost::filesystem::directory_iterator(directory)) {
        cslibs_ndt::Indexed<T, Size, Dim> indexed
                = YAML::LoadFile(entry.path().string()).as<cslibs_ndt::Indexed<T, Size, Dim>>();
        storage->insert(indexed.index_, indexed.data_);
    }

    return true;
}
}

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
