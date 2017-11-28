#ifndef CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
#define CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math_2d/linear/point.hpp>

#include <cslibs_utility/synchronized/wrap_around.hpp>

namespace cslibs_ndt {
template<std::size_t Dim>
class DistributionContainer;

template<std::size_t Dim>
class DistributionContainer {
public:
    using distribution_container_t = DistributionContainer<Dim>;
    using distribution_t = cslibs_math::statistics::Distribution<Dim>;
    using mutex_t        = std::mutex;
    using lock_t         = std::unique_lock<mutex_t>;
    enum Action {NONE, ALLOCATED, TOUCHED};

    using handle_t = cslibs_utility::synchronized::WrapAround<distribution_container_t>;
    using const_handle_t = cslibs_utility::synchronized::WrapAround<const distribution_container_t>;

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

    inline handle_t getHandle()
    {
        return handle_t(this, &data_mutex_);
    }

    inline const_handle_t getHandle() const
    {
        return const_handle_t(this, &data_mutex_);
    }

private:
    mutable mutex_t data_mutex_;
    Action          action_;
    distribution_t  data_;

};
}
#endif // CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
