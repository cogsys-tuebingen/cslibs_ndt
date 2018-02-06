#ifndef CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
#define CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_utility/synchronized/wrap_around.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>


namespace cslibs_ndt {
template<std::size_t Dim>
class Distribution {
public:
    /// insert a probabilistic decay function and model map update via bresenham
    /// use the standard occupancy grid map update model

    using distribution_container_t = Distribution<Dim>;
    using distribution_t = cslibs_math::statistics::Distribution<Dim, 3>;
    using mutex_t        = std::mutex;
    using lock_t         = std::unique_lock<mutex_t>;

    using handle_t       = cslibs_utility::synchronized::WrapAround<distribution_container_t>;
    using const_handle_t = cslibs_utility::synchronized::WrapAround<const distribution_container_t>;

    inline Distribution()
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

    inline void merge(const Distribution &)
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

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }

private:
    mutable mutex_t data_mutex_;
    distribution_t  data_;

};
}
#endif // CSLIBS_NDT_DISTRIBUTION_CONTAINER_HPP
