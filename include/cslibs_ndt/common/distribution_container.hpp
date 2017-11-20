#ifndef CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
#define CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>

namespace cslibs_ndt {
template<std::size_t Dim, bool limit_covariance = false>
class Distribution {
public:
    using distribution_t = cslibs_math::statistics::Distribution<Dim, limit_covariance>;
    using mutex_t        = std::mutex;

    Distribution() = default;
    virtual ~Distribution() = default;

    Distribution(const Distribution &other) :
        data_(other.data_)
    {
    }

    Distribution(Distribution &&other) :
        data_(std::move(other.data_))
    {
    }

    Distribution& operator = (const Distribution &other)
    {
        data_ = other.data_;
        return *this;
    }

    Distribution& operator = (Distribution &&other)
    {
        data_ = std::move(other.data_);
        return *this;
    }


    operator const vector_t& () const
    {
        return data_;
    }

    operator vector_t& ()
    {
        return data_;
    }

    operator vector_t () const
    {
        return data_;
    }

    operator vector_t* ()
    {
        return &data_;
    }

    inline const vector_t& data() const
    {
        return data_;
    }

    inline vector_t& data()
    {
        return data_;
    }


    inline void merge(const Chunk &)
    {
    }

    inline void lock() const
    {
        data_mutex_.lock();
    }

    inline void unlock() const
    {
        data_mutex_.unlock();
    }


private:
    mutex_t        data_mutex_;
    distribution_t data_;

};
}
#endif // CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
