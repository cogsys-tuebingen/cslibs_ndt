#ifndef CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math/statistics/stable_distribution.hpp>
#include <cslibs_gridmaps/utility/inverse_model.hpp>

#include <cslibs_indexed_storage/storage.hpp>

namespace cslibs_ndt {
template<typename T, std::size_t Dim>
class /*EIGEN_ALIGN16*/ OccupancyDistribution
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    using allocator_t               = Eigen::aligned_allocator<OccupancyDistribution<T,Dim>>;

    using Ptr                       = std::shared_ptr<OccupancyDistribution<T,Dim>>;
    using distribution_container_t  = OccupancyDistribution<T, Dim>;
    using distribution_t            = cslibs_math::statistics::StableDistribution<T,Dim,3>;
    using distribution_ptr_t        = typename distribution_t::Ptr;
    using point_t                   = typename distribution_t::sample_t;
    using ivm_t                     = cslibs_gridmaps::utility::InverseModel<T>;

    inline OccupancyDistribution() :
        num_free_(0ul)
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
        distribution_(other.distribution_ ? new distribution_t(*(other.distribution_)) : nullptr)
    {
    }

    inline OccupancyDistribution& operator = (const OccupancyDistribution &other)
    {
        num_free_ = other.num_free_;
        if (other.distribution_)
            distribution_.reset(new distribution_t(*(other.distribution_)));
        else
            distribution_.reset();
        return *this;
    }

    inline void updateFree()
    {
        ++num_free_;
    }

    inline void updateFree(const std::size_t &n)
    {
        num_free_ += n;
    }

    inline void updateOccupied(const point_t &p)
    {
        if (!distribution_)
            distribution_.reset(new distribution_t());

        *distribution_ += p;
    }

    inline void updateOccupied(const distribution_t &d)
    {
        if (!distribution_)
            distribution_.reset(new distribution_t(d));
        else
            *distribution_ += d;
    }

    inline std::size_t numFree() const
    {
        return num_free_;
    }

    inline std::size_t numOccupied() const
    {
        return distribution_ ? distribution_->getN() : 0ul;
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
                    num_free_ * inverse_model.getLogOddsFree() +
                    distribution_->getN() * inverse_model.getLogOddsOccupied() -
                    (num_free_ + distribution_->getN() - 1) * inverse_model.getLogOddsPrior())
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

    inline void merge(const OccupancyDistribution &other)
    {
        num_free_ += other.num_free_;
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
    std::size_t        num_free_     = 0ul;
    distribution_ptr_t distribution_ = nullptr;
};
}

#endif // CSLIBS_NDT_COMMON_OCCUPANCY_DISTRIBUTION_HPP
