#ifndef CSLIBS_NDT_3D_STATIC_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_STATIC_OCCUPANCY_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_3d/linear/pose.hpp>
#include <cslibs_math_3d/linear/point.hpp>

#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

#include <cslibs_math_3d/algorithms/bresenham.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace static_maps {
class OccupancyGridmap
{
public:
    using Ptr                               = std::shared_ptr<OccupancyGridmap>;
    using pose_t                            = cslibs_math_3d::Pose3d;
    using transform_t                       = cslibs_math_3d::Transform3d;
    using point_t                           = cslibs_math_3d::Point3d;
    using index_t                           = std::array<int, 3>;
    using size_t                            = std::array<std::size_t, 3>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::OccupancyDistribution<3>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 8>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 8>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 8>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::array::Array>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;
    using line_iterator_t                   = cslibs_math_3d::algorithms::Bresenham;

    OccupancyGridmap(const pose_t &origin,
            const double           resolution,
            const size_t          &size) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        bundle_resolution_2_(0.25 * bundle_resolution_ * bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0],
                                                                 size[1],
                                                                 size[2]);
        for(std::size_t i = 1 ; i < 8 ; ++i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1,
                                                                     size[1] + 1,
                                                                     size[2] + 1);
        }
        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2,
                                                                     size[1] * 2,
                                                                     size[2] * 2);
    }

    OccupancyGridmap(const double  origin_x,
                     const double  origin_y,
                     const double  origin_phi,
                     const double  resolution,
                     const size_t &size) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        bundle_resolution_2_(0.25 * bundle_resolution_ * bundle_resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t),
                  distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0],
                                                                 size[1],
                                                                 size[2]);
        for(std::size_t i = 1 ; i < 8 ; ++i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1,
                                                                     size[1] + 1,
                                                                     size[2] + 1);
        }
        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2,
                                                                     size[1] * 2,
                                                                     size[2] * 2);
    }

    inline pose_t getOrigin() const
    {
        return w_T_m_;
    }

    inline void add(const point_t &start_p,
                    const point_t &end_p)
    {
        const index_t start_index = toBundleIndex(start_p);
        const index_t end_index   = toBundleIndex(end_p);
        line_iterator_t it(start_index, end_index);

        while (!it.done()) {
            const index_t bi = {{it.x(), it.y(), it.z()}};
            (it.distance2() > bundle_resolution_2_) ?
                        updateFree(bi) :
                        updateOccupied(bi, end_p);
            ++ it;
        }
        updateOccupied(end_index, end_p);
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return getAllocate(bi);
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return getAllocate(bi);
    }

    inline double getBundleResolution() const
    {
        return bundle_resolution_;
    }

    inline double getResolution() const
    {
        return resolution_;
    }

    inline double getHeight() const
    {
        return size_[1] * resolution_;
    }

    inline double getWidth() const
    {
        return size_[0] * resolution_;
    }

    inline size_t getSize() const
    {
        return size_;
    }

    inline size_t getBundleSize() const
    {
        return {{size_[0] * 2, size_[1] * 2}};
    }

protected:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const double                                    bundle_resolution_;
    const double                                    bundle_resolution_inv_;
    const double                                    bundle_resolution_2_;
    const transform_t                               w_T_m_;
    const transform_t                               m_T_w_;
    const size_t                                    size_;

    mutable mutex_t                                 storage_mutex_;
    mutable distribution_storage_array_t            storage_;
    mutable mutex_t                                 bundle_storage_mutex_;
    mutable distribution_bundle_storage_ptr_t       bundle_storage_;

    inline distribution_t* getAllocate(const distribution_storage_ptr_t &s,
                                       const index_t &i) const
    {
        lock_t l(storage_mutex_);
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }

        auto allocate_bundle = [this, &bi]() {
            distribution_bundle_t b;
            const int divx = cslibs_math::common::div<int>(bi[0], 2);
            const int divy = cslibs_math::common::div<int>(bi[1], 2);
            const int divz = cslibs_math::common::div<int>(bi[2], 2);
            const int modx = cslibs_math::common::mod<int>(bi[0], 2);
            const int mody = cslibs_math::common::mod<int>(bi[1], 2);
            const int modz = cslibs_math::common::mod<int>(bi[2], 2);

            const index_t storage_0_index = {{divx,        divy,        divz}};
            const index_t storage_1_index = {{divx + modx, divy,        divz}};
            const index_t storage_2_index = {{divx,        divy + mody, divz}};
            const index_t storage_3_index = {{divx + modx, divy + mody, divz}};
            const index_t storage_4_index = {{divx,        divy,        divz + modz}};
            const index_t storage_5_index = {{divx + modx, divy,        divz + modz}};
            const index_t storage_6_index = {{divx,        divy + mody, divz + modz}};
            const index_t storage_7_index = {{divx + modx, divy + mody, divz + modz}};

            b[0] = getAllocate(storage_[0], storage_0_index);
            b[1] = getAllocate(storage_[1], storage_1_index);
            b[2] = getAllocate(storage_[2], storage_2_index);
            b[3] = getAllocate(storage_[3], storage_3_index);
            b[4] = getAllocate(storage_[4], storage_4_index);
            b[5] = getAllocate(storage_[5], storage_5_index);
            b[6] = getAllocate(storage_[6], storage_6_index);
            b[7] = getAllocate(storage_[7], storage_7_index);

            lock_t(bundle_storage_mutex_);
            return &(bundle_storage_->insert(bi, b));
        };

        return bundle == nullptr ? allocate_bundle() : bundle;
    }

    inline void updateFree(const index_t &bi) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->updateFree();
        bundle->at(1)->updateFree();
        bundle->at(2)->updateFree();
        bundle->at(3)->updateFree();
        bundle->at(4)->updateFree();
        bundle->at(5)->updateFree();
        bundle->at(6)->updateFree();
        bundle->at(7)->updateFree();
    }

    inline void updateOccupied(const index_t &bi,
                               const point_t &p) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->updateOccupied(p);
        bundle->at(1)->updateOccupied(p);
        bundle->at(2)->updateOccupied(p);
        bundle->at(3)->updateOccupied(p);
        bundle->at(4)->updateOccupied(p);
        bundle->at(5)->updateOccupied(p);
        bundle->at(6)->updateOccupied(p);
        bundle->at(7)->updateOccupied(p);
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
    }
};
}
}
#endif // CSLIBS_NDT_3D_STATIC_OCCUPANCY_GRIDMAP_HPP
