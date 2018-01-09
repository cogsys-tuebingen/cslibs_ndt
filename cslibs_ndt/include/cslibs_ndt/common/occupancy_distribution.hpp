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

    inline OccupancyDistribution(const std::size_t num_free,
                                 const std::size_t num_occ) :
        num_free_(num_free),
        num_occupied_(num_occ)
    {
    }

    inline OccupancyDistribution(const std::size_t    num_free,
                                 const std::size_t    num_occ,
                                 const distribution_t data) :
        num_free_(num_free),
        num_occupied_(num_occ),
        distribution_(new distribution_t(data))
    {
    }

    inline OccupancyDistribution(const OccupancyDistribution &other) :
        num_free_(other.num_free_),
        num_occupied_(other.num_occupied_),
      distribution_(other.distribution_)
    {
    }

    inline OccupancyDistribution& operator = (const OccupancyDistribution &other)
    {
        num_free_ = other.num_free_;
        num_occupied_ = other.num_occupied_;
        distribution_ = other.distribution_;
        return *this;
    }

    inline void updateFree()
    {
        ++ num_free_;
    }

    inline void updateOccupied(const point_t & p)
    {
        ++ num_occupied_;
        if (!distribution_)
            distribution_.reset(new distribution_t());

        distribution_->add(p);
    }

    inline std::size_t numFree() const
    {
        return num_free_;
    }

    inline std::size_t numOccupied() const
    {
        return num_occupied_;
    }

    inline const std::shared_ptr<distribution_t> getDistribution() const
    {
        return distribution_;
    }

    inline std::shared_ptr<distribution_t> getDistribution()
    {
        return distribution_;
    }

    inline void merge(const OccupancyDistribution&)
    {
    }

private:
    std::size_t                     num_free_;
    std::size_t                     num_occupied_;
    std::shared_ptr<distribution_t> distribution_;
} __attribute__ ((aligned (64)));
}

#endif // CSLIBS_NDT_OCCUPANCY_DISTRIBUTION_HPP
