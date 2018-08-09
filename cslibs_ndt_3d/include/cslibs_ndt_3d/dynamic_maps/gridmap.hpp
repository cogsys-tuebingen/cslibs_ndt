#ifndef CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP

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
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace dynamic_maps {
class Gridmap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Gridmap>;

    using Ptr                               = std::shared_ptr<Gridmap>;
    using pose_2d_t                         = cslibs_math_2d::Pose2d;
    using pose_t                            = cslibs_math_3d::Pose3d;
    using transform_t                       = cslibs_math_3d::Transform3d;
    using point_t                           = cslibs_math_3d::Point3d;
    using index_t                           = std::array<int, 3>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::Distribution<3>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 8>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 8>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 8>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;

    Gridmap(const pose_t        &origin,
            const double         resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
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
    }

    Gridmap(const pose_t &origin,
            const double &resolution,
            const index_t &min_index,
            const index_t &max_index,
            const std::shared_ptr<distribution_bundle_storage_t> &bundles,
            const distribution_storage_array_t                   &storage) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_(min_index),
        max_index_(max_index),
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline bool empty() const
    {
        return min_index_[0] == std::numeric_limits<int>::max();
    }

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    inline point_t getMin() const
    {
        lock_t(bundle_storage_mutex_);
        return point_t(min_index_[0] * bundle_resolution_,
                       min_index_[1] * bundle_resolution_,
                       min_index_[2] * bundle_resolution_);
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    inline point_t getMax() const
    {
        lock_t(bundle_storage_mutex_);
        return point_t((max_index_[0] + 1) * bundle_resolution_,
                       (max_index_[1] + 1) * bundle_resolution_,
                       (max_index_[2] + 1) * bundle_resolution_);
    }

    /**
     * @brief Get the origin of the map.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
        origin.translation() = getMin();
        return origin;
    }

    /**
     * @brief Get the origin of the map.
     * @return the initial origin
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

    inline void insert(const point_t &p,
                    index_t &bi)
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bi = toBundleIndex(p);
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
        for (const auto &p : *points) {
            const point_t pm = points_origin * p;
            if (pm.isNormal()) {
                const index_t &bi = toBundleIndex(pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->data().add(pm);
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


    inline index_t getMinBundleIndex() const
    {
        lock_t(bundle_storage_mutex_);
        return min_index_;
    }

    inline index_t getMaxDistributionIndex() const
    {
        lock_t(bundle_storage_mutex_);
        return max_index_;
    }

    inline const distribution_bundle_t* getDistributionBundle(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        return bundle;
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
        lock_t(bundle_storage_mutex_);
        return (max_index_[1] - min_index_[1] + 1) * bundle_resolution_;
    }

    inline double getWidth() const
    {
        lock_t(bundle_storage_mutex_);
        return (max_index_[0] - min_index_[0] + 1) * bundle_resolution_;
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
      lock_t l(bundle_storage_mutex_);
      const point_t p_m = m_T_w_ * p_w.translation();
      index_t i = {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_)),
                    static_cast<int>(std::floor(p_m(1) * bundle_resolution_)),
                    static_cast<int>(std::floor(p_m(2) * bundle_resolution_))}};

      return (i[0] >= min_index_[0]  && i[0] <= max_index_[0]) &&
             (i[1] >= min_index_[1]  && i[1] <= max_index_[1]) &&
             (i[2] >= min_index_[2]  && i[2] <= max_index_[2]);
    }

    inline virtual bool validate(const pose_2d_t &p_w) const
    {
      lock_t l(bundle_storage_mutex_);
      const point_t p_m = m_T_w_ * point_t(p_w.translation()(0), p_w.translation()(1), 0.0);
      index_t i = {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_)),
                    static_cast<int>(std::floor(p_m(1) * bundle_resolution_)),
                    0}};

      return (i[0] >= min_index_[0]  && i[0] <= max_index_[0]) &&
             (i[1] >= min_index_[1]  && i[1] <= max_index_[1]);
    }

protected:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const double                                    bundle_resolution_;
    const double                                    bundle_resolution_inv_;
    const transform_t                               w_T_m_;
    const transform_t                               m_T_w_;

    mutable index_t                                 min_index_;
    mutable index_t                                 max_index_;
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

    inline distribution_bundle_t* getAllocate(const index_t &bi) const
    {
        auto get_allocate = [this](const index_t &bi) {
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
                updateIndices(bi);
                return &(bundle_storage_->insert(bi, b));
            };
            return bundle ? bundle : allocate_bundle();
        };
        return get_allocate(bi);
    }

    inline void updateIndices(const index_t &chunk_index) const
    {
        min_index_ = std::min(min_index_, chunk_index);
        max_index_ = std::max(max_index_, chunk_index);
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;        
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
    }
}__attribute__ ((aligned (16)));
}
}

#endif // CSLIBS_NDT_3D_DYNAMIC_MAPS_GRIDMAP_HPP
