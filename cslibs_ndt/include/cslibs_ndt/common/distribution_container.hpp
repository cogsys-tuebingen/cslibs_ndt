#ifndef CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
#define CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP

#include <mutex>
#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math_2d/linear/point.hpp>

namespace cslibs_ndt {
template<std::size_t Dim>
class DistributionContainer {
public:
    using distribution_container_t = DistributionContainer<Dim>;
    using distribution_t = cslibs_math::statistics::Distribution<Dim>;
    using mutex_t        = std::mutex;
    using lock_t         = std::unique_lock<mutex_t>;
    enum Action {NONE, ALLOCATED, TOUCHED};

    class Handle {
    public:
        inline Handle() :
            distribution_(nullptr)
        {
        }

        inline Handle(distribution_container_t *container) :
            distribution_(container)
        {
            if(distribution_)
                distribution_->data_mutex_.lock();
        }

        virtual inline ~Handle()
        {
            if(distribution_)
               distribution_->data_mutex_.unlock();
        }

        inline bool empty() const
        {
            return distribution_ == nullptr;
        }

        inline const distribution_t& data() const
        {
            return distribution_->data_;
        }

        inline distribution_t& data()
        {
            return distribution_->data_;
        }

        inline operator distribution_t* ()
        {
            return distribution_;
        }

        inline operator distribution_t const *() const
        {
            return distribution_;
        }

        inline distribution_container_t * operator -> ()
        {
            return distribution_;
        }

        inline distribution_container_t const * operator -> () const
        {
            return distribution_;
        }

        inline operator distribution_container_t* ()
        {
            return distribution_;
        }

        inline operator distribution_container_t const *() const
        {
            return distribution_;
        }

    private:
        distribution_container_t   *distribution_;
    };

    using handle_t = Handle;
    friend class Handle;

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

private:
    Action          action_;
    mutable mutex_t data_mutex_;
    distribution_t  data_;

};
}
#endif // CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
