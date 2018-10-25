#ifndef CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<std::size_t Dim>
class EIGEN_ALIGN16 Distribution
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using allocator_t              = Eigen::aligned_allocator<Distribution>;
    using distribution_container_t = Distribution<Dim>;
    using distribution_t           = cslibs_math::statistics::Distribution<Dim, 3>;

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

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }

private:
    distribution_t  data_;
};
}
#endif // CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
