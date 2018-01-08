#ifndef CSLIBS_NDT_OCCUPANCY_DISTRIBUTION_HPP
#define CSLIBS_NDT_OCCUPANCY_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_utility/synchronized/wrap_around.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<std::size_t Dim>
class OccupancyDistribution {
public:
    using distribution_container_t  = OccupancyDistribution<Dim>;
    using distribution_t            = cslibs_math::statistics::Distribution<Dim, 3>;
    using point_t                   = typename distribution_t::sample_t;
    using mutex_t                   = std::mutex;
    using lock_t                    = std::unique_lock<mutex_t>;

    using handle_t                  = cslibs_utility::synchronized::WrapAround<distribution_container_t>;
    using const_handle_t            = cslibs_utility::synchronized::WrapAround<const distribution_container_t>;

    inline OccupancyDistribution() :
        num_free_(0),
        num_occupied_(0)
    {
    }

    inline OccupancyDistribution(const OccupancyDistribution &other) :
        data_(other.data_),
        num_free_(other.num_free_),
        num_occupied_(other.num_occupied_)
    {
    }

   inline OccupancyDistribution(OccupancyDistribution &&other) :
        data_(std::move(other.data_)),
        num_free_(std::move(other.num_free_)),
        num_occupied_(std::move(other.num_occupied_))
    {
    }

    inline ~OccupancyDistribution()
    {
        if (data_)
            delete data_;
    }

    inline OccupancyDistribution& operator = (const OccupancyDistribution &other)
    {
        data_ = other.data_;
        num_free_ = other.num_free_;
        num_occupied_ = other.num_occupied_;
        return *this;
    }

    inline OccupancyDistribution& operator = (OccupancyDistribution &&other)
    {
        data_ = std::move(other.data_);
        num_free_ = std::move(other.num_free_);
        num_occupied_ = std::move(other.num_occupied_);
        return *this;
    }

    inline void updateFree()
    {
        ++ num_free_;
    }

    inline void updateOccupied(const point_t & p)
    {
        ++ num_occupied_;
        if (!data_)
            data_ = new distribution_t();

        getHandle()->data()->add(p);
    }

    inline std::size_t numFree() const
    {
        return num_free_;
    }

    inline std::size_t numOccupied() const
    {
        return num_occupied_;
    }
/*
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
    }*/

    inline const distribution_t* data() const
    {
        return data_;
    }

    inline distribution_t* data()
    {
        return data_;
    }

    inline void merge(const OccupancyDistribution&)
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
    distribution_t* data_;

    std::size_t     num_free_;
    std::size_t     num_occupied_;
};
}

#endif // CSLIBS_NDT_OCCUPANCY_DISTRIBUTION_HPP
