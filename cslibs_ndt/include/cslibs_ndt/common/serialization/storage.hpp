#ifndef CSLIBS_NDT_STORAGE_HPP
#define CSLIBS_NDT_STORAGE_HPP

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_ndt/common/serialization/filesystem.hpp>
#include <cslibs_ndt/common/serialization/indexed_occupancy_distribution.hpp>
#include <cslibs_ndt/common/serialization/indexed_distribution.hpp>

#include <cslibs_ndt/common/indexed.hpp>
#include <cslibs_ndt/common/distribution.hpp>

#include <fstream>
#include <yaml-cpp/yaml.h>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {

template <template <std::size_t> class T, std::size_t Size>
void write(const T<Size> &d, std::ofstream &out)
{
    cslibs_math::serialization::distribution::binary<Size, 3>::write(d.data(), out);
}

template <template <std::size_t> class T, std::size_t Size>
std::size_t read(std::ifstream &in, T<Size> &d)
{
    return cslibs_math::serialization::distribution::binary<Size, 3>::read(in, d.data());
}

template<std::size_t Size>
void write(const OccupancyDistribution<Size> &d, std::ofstream &out)
{
    cslibs_math::serialization::io<std::size_t>::write(d.numFree(), out);
    if(!d.getDistribution()) {
        cslibs_math::serialization::distribution::binary<Size, 3>::write(out);
    } else{
        cslibs_math::serialization::distribution::binary<Size, 3>::write(*(d.getDistribution()), out);
    }

}

template<std::size_t Size>
std::size_t read(std::ifstream &in, OccupancyDistribution<Size> &d)
{
    std::size_t f = cslibs_math::serialization::io<std::size_t>::read(in);
    d = OccupancyDistribution<Size>(f);
    typename OccupancyDistribution<Size>::distribution_t tmp;
    std::size_t r = cslibs_math::serialization::distribution::binary<Size, 3>::read(in,tmp);
    if(tmp.getN() != 0) {
        d.getDistribution().reset(new typename OccupancyDistribution<Size>::distribution_t(tmp));
    }
    return  sizeof(std::size_t) + r;
}


template <template <std::size_t> class T, std::size_t Size, std::size_t Dim>
struct binary {
    using index_t   = std::array<int, Dim>;
    using data_t    = T<Size>;
    using storage_t =  cis::Storage<data_t, index_t, cis::backend::kdtree::KDTree>;

    inline static bool save(const std::shared_ptr<storage_t>   &storage,
                            const boost::filesystem::path      &path)
    {
        std::ofstream out(path.string(), std::ios::binary | std::ios::trunc);
        if(!out.is_open()) {
            std::cerr << "Could not open '" << path.string() << "'\n";
            return false;
        }

        auto write = [&out] (const index_t &index, const data_t &data) {
            cslibs_math::serialization::array::binary<int, Dim>::write(index, out);
            cslibs_ndt::write(data, out);
        };
        storage->traverse(write);
        out.close();
        return true;
    }

    inline static bool load(const boost::filesystem::path &path,
                            std::shared_ptr<storage_t> &storage)
    {
        storage.reset(new storage_t);
        std::ifstream in(path.string(), std::ios::binary);
        if(!in.is_open()) {
            std::cerr << "Could not open '" << path.string() << "'\n";
            return false;
        }

        try {
            in.seekg (0, std::ios::end);
            const std::size_t size = in.tellg();
            in.seekg (0, std::ios::beg);
            std::size_t read = 0;
            while(read < size) {
                index_t index;
                data_t  data;
                read += cslibs_math::serialization::array::binary<int, Dim>::read(in, index);
                read += cslibs_ndt::read(in, data);
                storage->insert(index, data);
            }
        } catch (const std::exception &e) {
            std::cerr << "Faild reading file '" << e.what() << "'\n";
            return false;
        }
        return true;
    }
};

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
