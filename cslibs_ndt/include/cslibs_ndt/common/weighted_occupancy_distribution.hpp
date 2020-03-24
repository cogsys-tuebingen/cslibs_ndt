#ifndef CSLIBS_NDT_COMMON_WEIGHTED_OCCUPANCY_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_WEIGHTED_OCCUPANCY_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/weighted_distribution.hpp>
#include <cslibs_math/statistics/stable_weighted_distribution.hpp>
#include <cslibs_gridmaps/utility/inverse_model.hpp>

#include <cslibs_indexed_storage/storage.hpp>

namespace cslibs_ndt {
template<typename T, std::size_t Dim>
class /*EIGEN_ALIGN16*/ WeightedOccupancyDistribution
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    using allocator_t               = Eigen::aligned_allocator<WeightedOccupancyDistribution<T,Dim>>;

    using Ptr                       = std::shared_ptr<WeightedOccupancyDistribution<T,Dim>>;
    using distribution_container_t  = WeightedOccupancyDistribution<T,Dim>;
    using distribution_t            = cslibs_math::statistics::StableWeightedDistribution<T,Dim,3>;
    using distribution_ptr_t        = typename distribution_t::Ptr;
    using point_t                   = typename distribution_t::sample_t;
    using ivm_t                     = cslibs_gridmaps::utility::InverseModel<T>;

    inline WeightedOccupancyDistribution() :
        weight_free_(0)
    {
    }

    inline WeightedOccupancyDistribution(const T weight_free) :
        weight_free_(weight_free)
    {
    }

    inline WeightedOccupancyDistribution(const T weight_free,
                                         const distribution_t data) :
        weight_free_(weight_free),
        distribution_(new distribution_t(data))
    {
    }

    inline WeightedOccupancyDistribution(const WeightedOccupancyDistribution &other) :
        weight_free_(other.weight_free_),
        distribution_(other.distribution_)
    {
    }

    inline WeightedOccupancyDistribution& operator = (const WeightedOccupancyDistribution &other)
    {
        weight_free_   = other.weight_free_;
        distribution_  = other.distribution_;
        return *this;
    }

    inline void updateFree(const T& weight_free = 1.0)
    {
        weight_free_  += weight_free;
    }

    inline void updateOccupied(const point_t& p, const T& w = 1.0)
    {
        if (!distribution_)
            distribution_.reset(new distribution_t());

        distribution_->add(p, w);
    }

    inline void updateOccupied(const distribution_ptr_t &d)
    {
        if (!d)
            return;

        if (!distribution_)
            distribution_.reset(new distribution_t());

        *distribution_ += *d;
    }

    inline T weightFree() const
    {
        return weight_free_;
    }

    inline T weightOccupied() const
    {
        return distribution_ ? distribution_->getWeight() : T();
    }

    inline T getOccupancy(const typename ivm_t::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("inverse model not set!");

        return getOccupancy(*inverse_model);
    }

    inline T getOccupancy(const ivm_t &inverse_model) const
    {
        return distribution_ ?
                    cslibs_math::common::LogOdds<T>::from(
                        weight_free_ * inverse_model.getLogOddsFree() +
                        distribution_->getWeight() * inverse_model.getLogOddsOccupied() -
                        (weight_free_ + distribution_->getWeight() - 1) * inverse_model.getLogOddsPrior())
                      : T(0.0);
    }

    inline const distribution_ptr_t &getDistribution() const
    {
        return distribution_;
    }

    inline distribution_ptr_t &getDistribution()
    {
        return distribution_;
    }

    inline void merge(const WeightedOccupancyDistribution &other)
    {
        weight_free_ += other.weight_free_;
        if (other.distribution_) {
            if (distribution_)
                *distribution_ += *(other.distribution_);
            else
                distribution_ = other.distribution_;
        }
    }

    inline std::size_t byte_size() const
    {
        return distribution_ ? (sizeof(*this) + sizeof(distribution_t)) : sizeof(*this);
    }

private:
    T                  weight_free_;
    distribution_ptr_t distribution_;
};
}

#endif // CSLIBS_NDT_COMMON_WEIGHTED_OCCUPANCY_DISTRIBUTION_HPP
