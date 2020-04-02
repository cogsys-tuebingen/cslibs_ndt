#ifndef CSLIBS_NDT_SERIALIZATION_STORAGE_HPP
#define CSLIBS_NDT_SERIALIZATION_STORAGE_HPP

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_ndt/common/weighted_occupancy_distribution.hpp>
#include <cslibs_ndt/serialization/filesystem.hpp>

#include <cslibs_math/serialization/array.hpp>
#include <cslibs_math/serialization/distribution.hpp>
#include <cslibs_math/serialization/weighted_distribution.hpp>
#include <cslibs_math/serialization/stable_distribution.hpp>
#include <cslibs_math/serialization/stable_weighted_distribution.hpp>

#include <fstream>
#include <yaml-cpp/yaml.h>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt {
template <typename Tp, std::size_t Size>
void write(const Distribution<Tp,Size> &d, std::ofstream &out)
{
    if (!d.getDistribution())
        cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::write(out);
    else
        cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::write(*(d.getDistribution()), out);
}

template <typename Tp, template <typename,std::size_t> class T, std::size_t Size>
std::size_t read(std::ifstream &in, T<Tp,Size> &d)
{
    typename Distribution<Tp,Size>::distribution_t tmp;
    std::size_t r = cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::read(in,tmp);
    if (tmp.getN() != 0)
        d.getDistribution().reset(new typename Distribution<Tp,Size>::distribution_t(tmp));
    return r;
}

template<typename Tp, std::size_t Size>
void write(const OccupancyDistribution<Tp,Size> &d, std::ofstream &out)
{
    cslibs_math::serialization::io<std::size_t>::write(d.numFree(), out);
    if (!d.getDistribution())
        cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::write(out);
    else
        cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::write(*(d.getDistribution()), out);
}

template<typename Tp, std::size_t Size>
std::size_t read(std::ifstream &in, OccupancyDistribution<Tp,Size> &d)
{
    std::size_t f = cslibs_math::serialization::io<std::size_t>::read(in);
    d = OccupancyDistribution<Tp,Size>(f);
    typename OccupancyDistribution<Tp,Size>::distribution_t tmp;
    std::size_t r = cslibs_math::serialization::binary<cslibs_math::statistics::StableDistribution,Tp,Size,3>::read(in,tmp);
    if (tmp.getN() != 0)
        d.getDistribution().reset(new typename OccupancyDistribution<Tp,Size>::distribution_t(tmp));
    return sizeof(std::size_t) + r;
}

template<typename Tp, std::size_t Size>
void write(const WeightedOccupancyDistribution<Tp,Size> &d, std::ofstream &out)
{
    cslibs_math::serialization::io<Tp>::write(d.weightFree(), out);
    if (!d.getDistribution())
        cslibs_math::serialization::binary<cslibs_math::statistics::StableWeightedDistribution,Tp,Size,3>::write(out);
    else
        cslibs_math::serialization::binary<cslibs_math::statistics::StableWeightedDistribution,Tp,Size,3>::write(*(d.getDistribution()), out);
}

template<typename Tp, std::size_t Size>
std::size_t read(std::ifstream &in, WeightedOccupancyDistribution<Tp,Size> &d)
{
    //std::size_t n = cslibs_math::serialization::io<std::size_t>::read(in);
    Tp f = cslibs_math::serialization::io<Tp>::read(in);
    using distr_t = typename WeightedOccupancyDistribution<Tp,Size>::distribution_t;
    distr_t tmp;
    std::size_t r = cslibs_math::serialization::binary<cslibs_math::statistics::StableWeightedDistribution,Tp,Size,3>::read(in,tmp);
    d = WeightedOccupancyDistribution<Tp,Size>(f);
    if (tmp.getSampleCount() > 0)
        d.getDistribution().reset(new distr_t(tmp));
    return sizeof(std::size_t) + sizeof(Tp) + r;
}

template <template <typename,std::size_t> class T, typename Tp, std::size_t Size, std::size_t Dim,
          template <typename, typename, typename...> class backend_t>
struct binary {
    using index_t   = std::array<int, Dim>;
    using size_t    = std::array<std::size_t, Dim>;
    using data_t    = T<Tp,Size>;
    using storage_t = cis::Storage<data_t, index_t, backend_t>;

    inline static bool save(const std::shared_ptr<storage_t> &storage,
                            const boost::filesystem::path &path)
    {
        std::ofstream out(path.string(), std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            std::cerr << "Could not open '" << path.string() << std::endl;
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
        return loadStorage(path, storage);
    }

    inline static bool load(const boost::filesystem::path &path,
                            std::shared_ptr<storage_t> &storage,
                            const size_t &size,
                            const index_t &offset)
    {
        storage.reset(new storage_t);
        storage->template set<cis::option::tags::array_size>(size);
        storage->template set<cis::option::tags::array_offset>(offset);
        return loadStorage(path, storage);
    }

private:
    inline static bool loadStorage(const boost::filesystem::path &path,
                                   std::shared_ptr<storage_t> &storage)
    {
        std::ifstream in(path.string(), std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "Could not open '" << path.string() << std::endl;
            return false;
        }

        try {
            in.seekg (0, std::ios::end);
            const std::size_t size = in.tellg();
            in.seekg (0, std::ios::beg);
            std::size_t read = 0;
            while (read < size) {
                index_t index;
                data_t  data;
                read += cslibs_math::serialization::array::binary<int, Dim>::read(in, index);
                read += cslibs_ndt::read(in, data);
                storage->insert(index, data);
            }
        } catch (const std::exception &e) {
            std::cerr << "Failed reading file '" << e.what() << std::endl;
            return false;
        }
        return true;
    }
};
}

#endif // CSLIBS_NDT_SERIALIZATION_STORAGE_HPP
