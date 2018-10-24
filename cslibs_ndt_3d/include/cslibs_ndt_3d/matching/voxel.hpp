#ifndef CSLIBS_NDT_3D_VOXEL_HPP
#define CSLIBS_NDT_3D_VOXEL_HPP

#include <cslibs_math_3d/linear/pointcloud.hpp>
#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace matching {
class EIGEN_ALIGN16 Voxel {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using index_t = std::array<int, 3>;
    using point_t = cslibs_math_3d::Point3d;

    inline Voxel() :
        n_(1),
        n_1_(0)
    {
    }

    inline Voxel(const point_t &pt) :
        n_(2),
        n_1_(1),
        mean_(pt)
    {
    }

    inline virtual ~Voxel() = default;

    inline Voxel(const Voxel &other) :
        mean_(other.mean_)
    {
    }

   inline  Voxel(Voxel &&other) :
        mean_(std::move(other.mean_))
    {
    }

    inline Voxel& operator = (const Voxel &other)
    {
        mean_ = other.mean_;
        return *this;
    }

    inline Voxel& operator = (Voxel &&other)
    {
        mean_ = std::move(other.mean_);
        return *this;
    }

    inline point_t const & mean() const
    {
        return mean_;
    }

    inline void merge(const Voxel &other)
    {
        const std::size_t   _n  = n_1_ + other.n_1_;
        const point_t       _pt = (mean_ * static_cast<double>(n_1_) + other.mean_ * static_cast<double>(other.n_1_)) / static_cast<double>(_n);
        n_                      = _n + 1;
        n_1_                    = _n;
        mean_                   = _pt;
    }

private:
    std::size_t n_;
    std::size_t n_1_;
    point_t     mean_;
};

inline Voxel::index_t getIndex(const Voxel::point_t &p, const double inverse_resolution)
{
    return {{static_cast<int>(std::floor(p(0) * inverse_resolution)),
             static_cast<int>(std::floor(p(1) * inverse_resolution)),
             static_cast<int>(std::floor(p(2) * inverse_resolution))}};
}

using StaticVoxelGrid  = cis::Storage<Voxel, Voxel::index_t, cis::backend::array::Array>;
}
}


#endif // CSLIBS_NDT_3D_VOXEL_HPP
