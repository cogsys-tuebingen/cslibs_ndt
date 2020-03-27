#ifndef CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math/statistics/stable_distribution.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<typename T, std::size_t Dim>
class /*EIGEN_ALIGN16*/ Distribution //: public cslibs_math::statistics::StableDistribution<T,Dim,3>
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    using allocator_t = Eigen::aligned_allocator<Distribution<T,Dim>>;

    using Ptr                = std::shared_ptr<Distribution<T,Dim>>;
    using distribution_t     = cslibs_math::statistics::StableDistribution<T,Dim,3>;
    using distribution_ptr_t = typename distribution_t::Ptr;

    inline Distribution() = default;

    inline Distribution(const distribution_t& data) :
        distribution_(new distribution_t(data))
    {
    }

    inline Distribution(const Distribution &other) :
        distribution_(other.distribution_ ? new distribution_t(*(other.distribution_)) : nullptr)
    {
    }

    inline  Distribution(Distribution &&other) :
        distribution_(other.distribution_ ? new distribution_t(*(other.distribution_)) : nullptr)
    {
    }

    inline Distribution& operator = (const Distribution &other)
    {
        if (other.distribution_)
            distribution_.reset(new distribution_t(*(other.distribution_)));
        else
            distribution_.reset();
        return *this;
    }

    inline Distribution& operator = (Distribution &&other)
    {
        if (other.distribution_)
            distribution_.reset(new distribution_t(*(other.distribution_)));
        else
            distribution_.reset();
        return *this;
    }

    inline std::size_t byte_size() const
    {
        return distribution_ ? (sizeof(*this) + sizeof(distribution_t)) : sizeof(*this);
    }

    inline void update(const distribution_t &d)
    {
        if (!distribution_)
            distribution_.reset(new distribution_t(d));
        else
            *distribution_ += d;
    }

    inline const distribution_ptr_t &getDistribution() const
    {
        return distribution_;
    }

    inline distribution_ptr_t &getDistribution()
    {
        return distribution_;
    }

    inline void merge(const Distribution &other)
    {
        if (other.distribution_) {
            if (distribution_)
                *distribution_ += *(other.distribution_);
            else
                distribution_ = other.distribution_;
        }
    }

private:
    distribution_ptr_t distribution_ = nullptr;
};
}

#endif // CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
