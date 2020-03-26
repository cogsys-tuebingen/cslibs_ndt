#ifndef CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
#define CSLIBS_NDT_COMMON_DISTRIBUTION_HPP

#include <mutex>

#include <cslibs_math/statistics/distribution.hpp>
#include <cslibs_math/statistics/stable_distribution.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cslibs_ndt {
template<typename T, std::size_t Dim>
class EIGEN_ALIGN16 Distribution : public cslibs_math::statistics::StableDistribution<T,Dim,3>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Distribution<T,Dim>>;
    using distribution_t = cslibs_math::statistics::StableDistribution<T,Dim,3>;

    inline const distribution_t& data() const
    {
        return *this;
    }

    inline distribution_t& data()
    {
        return *this;
    }

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }
};
}

#endif // CSLIBS_NDT_COMMON_DISTRIBUTION_HPP
