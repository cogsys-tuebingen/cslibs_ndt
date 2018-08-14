#ifndef CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_2d/linear/pose.hpp>

#include <cslibs_math_3d/linear/pose.hpp>
#include <cslibs_math_3d/linear/point.hpp>

#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/linear/pointcloud.hpp>
#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace static_maps {
class Gridmap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Gridmap>;

    using Ptr                               = std::shared_ptr<Gridmap>;
    using ConstPtr                          = std::shared_ptr<Gridmap>;
    using pose_2d_t                         = cslibs_math_2d::Pose2d;
    using pose_t                            = cslibs_math_3d::Pose3d;
    using transform_t                       = cslibs_math_3d::Transform3d;
    using point_t                           = cslibs_math_3d::Point3d;
    using index_t                           = std::array<int, 3>;
    using size_t                            = std::array<std::size_t, 3>;
    using size_m_t                          = std::array<double, 3>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::Distribution<3>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 8>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 8>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 8>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::array::Array>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;

    Gridmap(const pose_t &origin,
            const double &resolution,
            const size_t &size,
            const index_t &min_bundle_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution,
                 (size[2] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
                           min_bundle_index[1] + static_cast<int>(size[1] * 2),
                           min_bundle_index[2] + static_cast<int>(size[2] * 2)}},
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
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1], size[2]);
        storage_[0]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                                                                   min_bundle_index[1] / 2,
                                                                   min_bundle_index[2] / 2);

        for(std::size_t i = 1 ; i < 8 ; ++ i) {
          storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1, size[2] + 1);
          storage_[i]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                                                                     min_bundle_index[1] / 2,
                                                                     min_bundle_index[2] / 2);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2, size[2] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0],
                                                                       min_bundle_index[1],
                                                                       min_bundle_index[2]);
    }

    Gridmap(const double &origin_x,
            const double &origin_y,
            const double &origin_phi,
            const double &resolution,
            const size_t &size,
            const index_t &min_bundle_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution,
                 (size[2] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
                           min_bundle_index[1] + static_cast<int>(size[1] * 2),
                           min_bundle_index[2] + static_cast<int>(size[2] * 2)}},
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
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1], size[2]);
        storage_[0]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                                                                   min_bundle_index[1] / 2,
                                                                   min_bundle_index[2] / 2);

        for(std::size_t i = 1 ; i < 8 ; ++ i) {
          storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1, size[2] + 1);
          storage_[i]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                                                                     min_bundle_index[1] / 2,
                                                                     min_bundle_index[2] / 2);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2, size[2] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0],
                                                                       min_bundle_index[1],
                                                                       min_bundle_index[2]);
    }

    Gridmap(const pose_t &origin,
            const double &resolution,
            const size_t &size,
            const std::shared_ptr<distribution_bundle_storage_t> &bundles,
            const distribution_storage_array_t                   &storage,
            const index_t &min_bundle_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution,
                 (size[2] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
                           min_bundle_index[1] + static_cast<int>(size[1] * 2),
                           min_bundle_index[2] + static_cast<int>(size[2] * 2)}},
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    inline point_t getMin() const
    {
      return point_t(min_bundle_index_[0] * bundle_resolution_,
                     min_bundle_index_[1] * bundle_resolution_,
                     min_bundle_index_[2] * bundle_resolution_);
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    inline point_t getMax() const
    {
      return point_t((max_bundle_index_[0] + 1) * bundle_resolution_,
                     (max_bundle_index_[1] + 1) * bundle_resolution_,
                     (max_bundle_index_[2] + 1) * bundle_resolution_);
    }

    /**
     * @brief Get the origin.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
        origin.translation() = getMin();
        return origin;
    }

    /**
     * @brief Get the initial origin of the map.
     * @return the inital origin
     */
    inline pose_t getInitialOrigin() const
    {
        return w_T_m_;
    }

    inline void insert(const point_t &p)
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            const index_t bi = toBundleIndex(p);
            if(!valid(bi))
                return;

            bundle = getAllocate(bi);
        }

        bundle->at(0)->getHandle()->data().add(p);
        bundle->at(1)->getHandle()->data().add(p);
        bundle->at(2)->getHandle()->data().add(p);
        bundle->at(3)->getHandle()->data().add(p);
        bundle->at(4)->getHandle()->data().add(p);
        bundle->at(5)->getHandle()->data().add(p);
        bundle->at(6)->getHandle()->data().add(p);
        bundle->at(7)->getHandle()->data().add(p);
    }

    inline void insert(const typename cslibs_math::linear::Pointcloud<point_t>::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        distribution_storage_t storage;
        storage.template set<cis::option::tags::array_size>(size_[0] * 2, size_[1] * 2, size_[2] * 2);
        storage.template set<cis::option::tags::array_offset>(min_bundle_index_[0],
                                                              min_bundle_index_[1],
                                                              min_bundle_index_[2]);
        for (const auto &p : *points) {
            const point_t pm = points_origin * p;
            if (pm.isNormal()) {
                const index_t &bi = toBundleIndex(pm);
                if(valid(bi)) {
                  distribution_t *d = storage.get(bi);
                  (d ? d : &storage.insert(bi, distribution_t()))->data().add(pm);
                }
            }
        }

        storage.traverse([this](const index_t& bi, const distribution_t &d) {
            distribution_bundle_t *bundle;
            {
                lock_t(bundle_storage_mutex_);
                bundle = getAllocate(bi);
            }
            bundle->at(0)->getHandle()->data() += d.data();
            bundle->at(1)->getHandle()->data() += d.data();
            bundle->at(2)->getHandle()->data() += d.data();
            bundle->at(3)->getHandle()->data() += d.data();
            bundle->at(4)->getHandle()->data() += d.data();
            bundle->at(5)->getHandle()->data() += d.data();
            bundle->at(6)->getHandle()->data() += d.data();
            bundle->at(7)->getHandle()->data() += d.data();
        });
    }

    inline double sample(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        if(!valid(bi))
            return 0.0;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto evaluate = [&p, &bundle]() {
            return 0.125 * (bundle->at(0)->getHandle()->data().sample(p) +
                            bundle->at(1)->getHandle()->data().sample(p) +
                            bundle->at(2)->getHandle()->data().sample(p) +
                            bundle->at(3)->getHandle()->data().sample(p) +
                            bundle->at(4)->getHandle()->data().sample(p) +
                            bundle->at(5)->getHandle()->data().sample(p) +
                            bundle->at(6)->getHandle()->data().sample(p) +
                            bundle->at(7)->getHandle()->data().sample(p));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline double sampleNonNormalized(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        if(!valid(bi))
            return 0.0;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto evaluate = [&p, &bundle]() {
            return 0.125 * (bundle->at(0)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(1)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(2)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(3)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(4)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(5)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(6)->getHandle()->data().sampleNonNormalized(p) +
                            bundle->at(7)->getHandle()->data().sampleNonNormalized(p));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        if(!valid(bi))
            return nullptr;

        lock_t l(bundle_storage_mutex_);
        return getAllocate(bi);
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        if(!valid(bi))
            return nullptr;

        lock_t l(bundle_storage_mutex_);
        return getAllocate(bi);
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        if(!valid(bi))
            return nullptr;

        distribution_bundle_t *bundle;
        {
            lock_t l(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        return bundle;
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
        return size_m_[1];
    }

    inline double getWidth() const
    {
        return size_m_[0];
    }

    inline size_t getSize() const
    {
        return size_;
    }

    inline size_t getBundleSize() const
    {
        return {{size_[0] * 2, size_[1] * 2, size_[2] * 2}};
    }

    inline distribution_storage_array_t const & getStorages() const
    {
        return storage_;
    }

    template <typename Fn>
    inline void traverse(const Fn& function) const
    {
        lock_t(storage_mutex_);
        lock_t(bundle_storage_mutex_);
        return bundle_storage_->traverse(function);
    }

    inline void getBundleIndices(std::vector<index_t> &indices) const
    {
        lock_t(storage_mutex_);
        lock_t(bundle_storage_mutex_);
        auto add_index = [&indices](const index_t &i, const distribution_bundle_t &d) {
            indices.emplace_back(i);
        };
        bundle_storage_->traverse(add_index);
    }

    inline std::size_t getByteSize() const
    {
        lock_t(storage_mutex_);
        lock_t(bundle_storage_mutex_);
        return sizeof(*this) +
                bundle_storage_->byte_size() +
                storage_[0]->byte_size() +
                storage_[1]->byte_size() +
                storage_[2]->byte_size() +
                storage_[3]->byte_size() +
                storage_[4]->byte_size() +
                storage_[5]->byte_size() +
                storage_[6]->byte_size() +
                storage_[7]->byte_size();
    }

    inline virtual bool validate(const pose_t &p_w) const
    {
      const point_t p_m = m_T_w_ * p_w.translation();
      const point_t min = getMin();
      const point_t max = getMax();
      return p_m(0) >= min(0) && p_m(0) < max(0) &&
             p_m(1) >= min(1) && p_m(1) < max(1) &&
             p_m(2) >= min(2) && p_m(2) < max(2);
    }

    inline virtual bool validate(const pose_2d_t &p_w) const
    {
      const point_t p_m = m_T_w_ * point_t(p_w.translation()(0), p_w.translation()(1), 0.0);
      const point_t min = getMin();
      const point_t max = getMax();
      return p_m(0) >= min(0) && p_m(0) < max(0) &&
             p_m(1) >= min(1) && p_m(1) < max(1);
    }

protected:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const double                                    bundle_resolution_;
    const double                                    bundle_resolution_inv_;
    const transform_t                               w_T_m_;
    const transform_t                               m_T_w_;
    const size_t                                    size_;
    const size_m_t                                  size_m_;
    const index_t                                   min_bundle_index_;
    const index_t                                   max_bundle_index_;

    mutable mutex_t                                 storage_mutex_;
    mutable distribution_storage_array_t            storage_;
    mutable mutex_t                                 bundle_storage_mutex_;
    mutable distribution_bundle_storage_ptr_t       bundle_storage_;

    inline distribution_t* getAllocate(const distribution_storage_ptr_t &s,
                                       const index_t &i) const
    {
        lock_t(storage_mutex_);
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
    {
        auto get_allocate = [this](const index_t &bi) {
            distribution_bundle_t *bundle = bundle_storage_->get(bi);

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
            return bundle ? bundle : allocate_bundle();
        };
        return get_allocate(bi);
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
    }

    inline bool valid(const index_t &bi) const
    {
        return bi[0] >= min_bundle_index_[0] &&
               bi[1] >= min_bundle_index_[1] &&
               bi[2] >= min_bundle_index_[2] &&
               bi[0] <= max_bundle_index_[0] &&
               bi[1] <= max_bundle_index_[1] &&
               bi[2] <= max_bundle_index_[2];

    }

}__attribute__ ((aligned (16)));
}
}

#endif // CSLIBS_NDT_3D_STATIC_MAPS_GRIDMAP_HPP
