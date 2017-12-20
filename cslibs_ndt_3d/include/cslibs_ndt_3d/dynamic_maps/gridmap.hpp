#ifndef CSLIBS_NDT_3D_DYNAMIC_GRIDMAP_HPP
#define CSLIBS_NDT_3D_DYNAMIC_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_3d/linear/pose.hpp>
#include <cslibs_math_3d/linear/point.hpp>

#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

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
    using Ptr                               = std::shared_ptr<Gridmap>;
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
    using offest_array_t                    = std::array<point_t, 7>;
    using index_array_t                     = std::array<index_t, 8>;
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
        offsets_{{point_t(-bundle_resolution_, 0.0, 0.0),
                  point_t(0.0, -bundle_resolution_, 0.0),
                  point_t(-bundle_resolution_, -bundle_resolution_, 0.0),
                  point_t(0.0, 0.0, -bundle_resolution_),
                  point_t(-bundle_resolution_, 0.0, -bundle_resolution_),
                  point_t(0.0, -bundle_resolution_, -bundle_resolution_),
                  point_t(-bundle_resolution_)}},
        index_offsets_{{index_t({0, 0, 0}),
                       index_t({-1, 0, 0}),
                       index_t({0, -1, 0}),
                       index_t({-1, -1, 0}),
                       index_t({0, 0, -1}),
                       index_t({-1, 0, -1}),
                       index_t({0, -1, -1}),
                       index_t({-1, -1, -1})}},
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

    inline point_t getMin() const
    {
        lock_t l(bundle_storage_mutex_);
        return point_t(min_index_[0] * bundle_resolution_,
                       min_index_[1] * bundle_resolution_,
                       min_index_[2] * bundle_resolution_);
    }

    inline point_t getMax() const
    {
        lock_t l(bundle_storage_mutex_);
        return point_t((max_index_[0] + 1) * bundle_resolution_,
                       (max_index_[1] + 1) * bundle_resolution_,
                       (max_index_[2] + 1) * bundle_resolution_);
    }

    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
        origin.translation() = getMin();
        return origin;
    }

    inline pose_t getInitialOrigin() const
    {
        return w_T_m_;
    }

    inline void add(const point_t &p)
    {
        distribution_bundle_t *bundle;
        {
            lock_t l(bundle_storage_mutex_);
            const index_t bi = toBundleIndex(p);
            bundle = bundle_storage_->get(bi);
            if(!bundle) {
                distribution_bundle_t b;
                b[0] = getAllocate(storage_[0], toIndex(p));
                b[1] = getAllocate(storage_[1], toIndex(p, offsets_[0]));
                b[2] = getAllocate(storage_[2], toIndex(p, offsets_[1]));
                b[3] = getAllocate(storage_[3], toIndex(p, offsets_[2]));
                b[4] = getAllocate(storage_[4], toIndex(p, offsets_[3]));
                b[5] = getAllocate(storage_[5], toIndex(p, offsets_[4]));
                b[6] = getAllocate(storage_[6], toIndex(p, offsets_[5]));
                b[7] = getAllocate(storage_[7], toIndex(p, offsets_[6]));
                bundle = &(bundle_storage_->insert(bi, b));
            }
            updateIndices(bi);
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

    inline void add(const point_t &p,
                    index_t &bi)
    {
        distribution_bundle_t *bundle;
        bi = toBundleIndex(p);
        {
            lock_t l(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
            if(!bundle) {
                distribution_bundle_t b;
                const index_array_t indices = toIndices(bi);
                b[0] = getAllocate(storage_[0], indices[0]);
                b[1] = getAllocate(storage_[1], indices[1]);
                b[2] = getAllocate(storage_[2], indices[2]);
                b[3] = getAllocate(storage_[3], indices[3]);
                b[4] = getAllocate(storage_[4], indices[4]);
                b[5] = getAllocate(storage_[5], indices[5]);
                b[6] = getAllocate(storage_[6], indices[6]);
                b[7] = getAllocate(storage_[7], indices[7]);
                bundle = &(bundle_storage_->insert(bi, b));
            }
            updateIndices(bi);
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

    inline double sample(const point_t &p) const
    {
        const index_t bi = toBundleIndex(p);
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto evaluate = [&p, &bundle]() {
            return 0.125 * (bundle->at(0)->data().sample(p) +
                            bundle->at(1)->data().sample(p) +
                            bundle->at(2)->data().sample(p) +
                            bundle->at(3)->data().sample(p) +
                            bundle->at(4)->data().sample(p) +
                            bundle->at(5)->data().sample(p) +
                            bundle->at(6)->data().sample(p) +
                            bundle->at(7)->data().sample(p));
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
            return 0.125 * (bundle->at(0)->data().sampleNonNormalized(p) +
                            bundle->at(1)->data().sampleNonNormalized(p) +
                            bundle->at(2)->data().sampleNonNormalized(p) +
                            bundle->at(3)->data().sampleNonNormalized(p) +
                            bundle->at(4)->data().sampleNonNormalized(p) +
                            bundle->at(5)->data().sampleNonNormalized(p) +
                            bundle->at(6)->data().sampleNonNormalized(p) +
                            bundle->at(7)->data().sampleNonNormalized(p));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline index_t getMinDistributionIndex() const
    {
        lock_t l(storage_mutex_);
        return min_index_;
    }

    inline index_t getMaxDistributionIndex() const
    {
        lock_t l(storage_mutex_);
        return max_index_;
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto allocate_bundle = [this, &bi]() {
            distribution_bundle_t b;
            const index_array_t indices = toIndices(bi);
            b[0] = getAllocate(storage_[0], indices[0]);
            b[1] = getAllocate(storage_[1], indices[1]);
            b[2] = getAllocate(storage_[2], indices[2]);
            b[3] = getAllocate(storage_[3], indices[3]);
            b[4] = getAllocate(storage_[4], indices[4]);
            b[5] = getAllocate(storage_[5], indices[5]);
            b[6] = getAllocate(storage_[6], indices[6]);
            b[7] = getAllocate(storage_[7], indices[7]);
            lock_t(bundle_storage_mutex_);
            updateIndices(bi);
            return &(bundle_storage_->insert(bi, b));
        };
        return bundle ? bundle : allocate_bundle();
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto allocate_bundle = [this, &bi]() {
            distribution_bundle_t b;
            const index_array_t indices = toIndices(bi);
            b[0] = getAllocate(storage_[0], indices[0]);
            b[1] = getAllocate(storage_[1], indices[1]);
            b[2] = getAllocate(storage_[2], indices[2]);
            b[3] = getAllocate(storage_[3], indices[3]);
            b[4] = getAllocate(storage_[4], indices[4]);
            b[5] = getAllocate(storage_[5], indices[5]);
            b[6] = getAllocate(storage_[6], indices[6]);
            b[7] = getAllocate(storage_[7], indices[7]);
            lock_t(bundle_storage_mutex_);
            updateIndices(bi);
            return &(bundle_storage_->insert(bi, b));
        };
        return bundle ? bundle : allocate_bundle();
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
        return (max_index_[1] - min_index_[1] + 1) * bundle_resolution_;
    }

    inline double getWidth() const
    {
        return (max_index_[0] - min_index_[0] + 1) * bundle_resolution_;
    }

protected:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const double                                    bundle_resolution_;
    const double                                    bundle_resolution_inv_;
    const transform_t                               w_T_m_;
    const transform_t                               m_T_w_;
    const offest_array_t                            offsets_;
    const index_array_t                             index_offsets_;

    mutable index_t                                 min_index_;
    mutable index_t                                 max_index_;
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

    inline void updateIndices(const index_t &chunk_index) const
    {
        min_index_ = std::min(min_index_, chunk_index);
        max_index_ = std::max(max_index_, chunk_index);
    }

    inline index_t toIndex(const point_t &p_w,
                           const point_t &off = point_t()) const
    {
        const point_t p_m = (m_T_w_ * p_w) + off;        
        return {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * resolution_inv_))}};
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;        
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
    }

    inline index_array_t toIndices(const index_t &bi) const
    {
        const index_t base_index({cslibs_math::common::div(bi[0], 2),
                                  cslibs_math::common::div(bi[1], 2),
                                  cslibs_math::common::div(bi[2], 2)});
        const index_t rel_index({1 - cslibs_math::common::mod(bi[0], 2),
                                 1 - cslibs_math::common::mod(bi[1], 2),
                                 1 - cslibs_math::common::mod(bi[2], 2)});
        return index_array_t({base_index + index_t({rel_index[0] * index_offsets_[0][0],
                                                    rel_index[1] * index_offsets_[0][1],
                                                    rel_index[2] * index_offsets_[0][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[1][0],
                                                    rel_index[1] * index_offsets_[1][1],
                                                    rel_index[2] * index_offsets_[1][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[2][0],
                                                    rel_index[1] * index_offsets_[2][1],
                                                    rel_index[2] * index_offsets_[2][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[3][0],
                                                    rel_index[1] * index_offsets_[3][1],
                                                    rel_index[2] * index_offsets_[3][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[4][0],
                                                    rel_index[1] * index_offsets_[4][1],
                                                    rel_index[2] * index_offsets_[4][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[5][0],
                                                    rel_index[1] * index_offsets_[5][1],
                                                    rel_index[2] * index_offsets_[5][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[6][0],
                                                    rel_index[1] * index_offsets_[6][1],
                                                    rel_index[2] * index_offsets_[6][2]}),
                              base_index + index_t({rel_index[0] * index_offsets_[7][0],
                                                    rel_index[1] * index_offsets_[7][1],
                                                    rel_index[2] * index_offsets_[7][2]})});
    }
};
}
}
#endif // CSLIBS_NDT_3D_DYNAMIC_GRIDMAP_HPP
