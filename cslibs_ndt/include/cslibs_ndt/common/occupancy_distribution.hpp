#ifndef CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_gridmaps/utility/inverse_model.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<std::size_t Dim>
class OccupancyDistribution
{
public:
    using Ptr                       = std::shared_ptr<OccupancyDistribution<Dim>>;
    using distribution_container_t  = OccupancyDistribution<Dim>;
    using distribution_t            = cslibs_math::statistics::Distribution<Dim, 3>;
    using distribution_ptr_t        = typename distribution_t::Ptr;
    using point_t                   = typename distribution_t::sample_t;
    using mutex_t                   = std::mutex;
    using lock_t                    = std::unique_lock<mutex_t>;

    inline OccupancyDistribution() :
        num_free_(0)
    {
    }

    inline OccupancyDistribution(const std::size_t num_free) :
        num_free_(num_free)
    {
    }

    inline OccupancyDistribution(const std::size_t    num_free,
                                 const distribution_t data) :
        num_free_(num_free),
        distribution_(new distribution_t(data))
    {
    }

    inline OccupancyDistribution(const OccupancyDistribution &other) :
        num_free_(other.num_free_),
        distribution_(other.distribution_),
        occupancy_(other.occupancy_),
        inverse_model_(other.inverse_model_)
    {
    }

    inline OccupancyDistribution& operator = (const OccupancyDistribution &other)
    {
        num_free_      = other.num_free_;
        distribution_  = other.distribution_;
        occupancy_     = other.occupancy_;
        inverse_model_ = other.inverse_model_;
        return *this;
    }

    inline void updateFree()
    {
        ++ num_free_;
        inverse_model_ = nullptr;
    }

    inline void updateFree(const std::size_t &num_free)
    {
        num_free_ += num_free;
        inverse_model_ = nullptr;
    }

    inline void updateOccupied(const point_t & p)
    {
        if (!distribution_)
            distribution_.reset(new distribution_t());

        distribution_->add(p);
        inverse_model_ = nullptr;
    }

    inline void updateOccupied(const distribution_ptr_t &d)
    {
        if (!d)
            return;

        if (!distribution_)
            distribution_.reset(new distribution_t());

        *distribution_ += *d;
        inverse_model_ = nullptr;
    }

    inline std::size_t numFree() const
    {
        return num_free_;
    }

    inline std::size_t numOccupied() const
    {
        return distribution_ ? distribution_->getN() : 0ul;
    }

    inline double getOccupancy(const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("inverse model not set!");

        if (inverse_model == inverse_model_)
            return occupancy_;

        inverse_model_ = inverse_model;
        occupancy_ = distribution_ ?
                    cslibs_math::common::LogOdds::from(
                        num_free_ * inverse_model_->getLogOddsFree() +
                        distribution_->getN() * inverse_model_->getLogOddsOccupied() -
                        (num_free_ + distribution_->getN()) * inverse_model_->getLogOddsPrior()) :
                    cslibs_math::common::LogOdds::from(
                        num_free_ * inverse_model_->getLogOddsFree() -
                        num_free_ * inverse_model_->getLogOddsPrior());
        return occupancy_;
    }

    inline const distribution_ptr_t &getDistribution() const
    {
        return distribution_;
    }

    inline distribution_ptr_t &getDistribution()
    {
        return distribution_;
    }

    inline void merge(const OccupancyDistribution&)
    {
    }

    inline std::size_t byte_size() const
    {
        return distribution_ ? (sizeof(*this) + sizeof(distribution_t)) : sizeof(*this);
    }

private:
    std::size_t        num_free_;
    distribution_ptr_t distribution_;

    mutable double                                      occupancy_;
    mutable cslibs_gridmaps::utility::InverseModel::Ptr inverse_model_;
} __attribute__ ((aligned (64)));
}

#endif // CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP
