#ifndef CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
#define CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP

#include <mutex>
#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math_2d/linear/point.hpp>

namespace cslibs_ndt {
template<std::size_t Dim, bool limit_covariance = false>
class DistributionContainer {
public:
    using distribution_t = cslibs_math::statistics::Distribution<Dim, limit_covariance>;
    using mutex_t        = std::mutex;
    using lock_t         = std::unique_lock<mutex_t>;

    enum Action {NONE, ALLOCATED, TOUCHED};

    inline DistributionContainer() :
        action_(ALLOCATED)
    {
    }

    inline virtual ~DistributionContainer() = default;

    inline DistributionContainer(const DistributionContainer &other) :
        data_(other.data_)
    {
    }

   inline  DistributionContainer(DistributionContainer &&other) :
        data_(std::move(other.data_))
    {
    }

    inline DistributionContainer& operator = (const DistributionContainer &other)
    {
        data_ = other.data_;
        return *this;
    }

    inline DistributionContainer& operator = (DistributionContainer &&other)
    {
        data_ = std::move(other.data_);
        return *this;
    }

    inline void setTouched()
    {
        if(action_ != ALLOCATED)
            action_ = TOUCHED;
    }

    inline void setNone()
    {
        action_ = NONE;
    }

    inline Action getAction() const
    {
        return action_;
    }

    inline operator const distribution_t& () const
    {
        return data_;
    }

    inline double sample(const cslibs_math_2d::Point2d &p) const
    {
        return data_.sample(p);
    }

    inline double sampleNonNormalized(const cslibs_math_2d::Point2d &p) const
    {
        return data_.sampleNonNormalized(p);
    }

    inline void add(const cslibs_math_2d::Point2d &p)
    {
        data_.add(p);
    }

    inline operator distribution_t& ()
    {
        return data_;
    }

    inline operator distribution_t () const
    {
        return data_;
    }

    inline operator distribution_t* ()
    {
        return &data_;
    }

    inline const distribution_t& data() const
    {
        return data_;
    }

    inline distribution_t& data()
    {
        return data_;
    }

    inline void merge(const DistributionContainer &)
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
    Action          action_;
    mutable mutex_t data_mutex_;
    distribution_t  data_;

};
}
#endif // CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
