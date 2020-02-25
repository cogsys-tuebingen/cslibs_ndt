#ifndef CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math/statistics/stable_distribution.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<typename T, std::size_t Dim>
class Distribution
{
public:
    using distribution_t = cslibs_math::statistics::StableDistribution<T,Dim,3>;

    inline Distribution()
    {
    }
    inline Distribution(const distribution_t& d) :
        data_(std::move(d))
    {
    }

    inline virtual ~Distribution() = default;

    inline Distribution(const Distribution &other) :
        data_(other.data_)
    {
    }

   inline  Distribution(Distribution &&other) :
        data_(std::move(other.data_))
    {
    }

    inline Distribution& operator = (const Distribution &other)
    {
        data_ = other.data_;
        return *this;
    }

    inline Distribution& operator = (Distribution &&other)
    {
        data_ = std::move(other.data_);
        return *this;
    }

    inline const distribution_t& data() const
    {
        return data_;
    }

    inline distribution_t& data()
    {
        return data_;
    }

    inline void merge(const Distribution &other)
    {
        data_ += other.data_;
    }

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }

private:
    distribution_t  data_;
};
}
#endif // CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
