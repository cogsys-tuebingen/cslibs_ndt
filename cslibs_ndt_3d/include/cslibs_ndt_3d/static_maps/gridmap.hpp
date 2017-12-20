#ifndef CSLIBS_NDT_3D_STATIC_GRIDMAP_HPP
#define CSLIBS_NDT_3D_STATIC_GRIDMAP_HPP

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
#include <cslibs_indexed_storage/backend/array/array.hpp>

#include <cslibs_utility/synchronized/wrap_around.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace static_maps {
class Gridmap
{
public:
    using Ptr                               = std::shared_ptr<Gridmap>;
    using pose_t                            = cslibs_math_3d::Pose3d;
    using transform_t                       = cslibs_math_3d::Transform3d;
    using point_t                           = cslibs_math_3d::Point3d;
    using index_t                           = std::array<int, 3>;
    using size_t                            = std::array<double, 3>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::Distribution<3>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 8>;
    using offest_array_t                    = std::array<point_t, 8>;
    using index_array_t                     = std::array<index_t, 8>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 8>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 8>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::array::Array>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;

    Gridmap(const pose_t        &origin,
            const double         resolution,
            const size_t        &size) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        offsets_{{point_t(),
                  point_t(-bundle_resolution_, 0.0, 0.0),
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
        min_index_{{static_cast<int>(std::floor(origin.tx() * bundle_resolution_inv_)),
                    static_cast<int>(std::floor(origin.ty() * bundle_resolution_inv_)),
                    static_cast<int>(std::floor(origin.tz() * bundle_resolution_inv_))}},
        max_index_{{static_cast<int>(std::floor((origin.tx() + size_[0]) * bundle_resolution_inv_)),
                    static_cast<int>(std::floor((origin.ty() + size_[1]) * bundle_resolution_inv_)),
                    static_cast<int>(std::floor((origin.tz() + size_[2]) * bundle_resolution_inv_))}},
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
        for(std::size_t i = 0 ; i < 8 ; ++i) {
            const point_t min    = w_T_m_.translation() + offsets_[i];
            const point_t max    = w_T_m_.translation() - offsets_[i];
            distribution_storage_ptr_t &storage = storage_[i];
            const index_t min_index = {{static_cast<int>(std::floor(min(0) * resolution_inv_)),
                                        static_cast<int>(std::floor(min(1) * resolution_inv_)),
                                        static_cast<int>(std::floor(min(2) * resolution_inv_))}};
            const index_t max_index = {{static_cast<int>(std::floor(max(0) * resolution_inv_)),
                                        static_cast<int>(std::floor(max(1) * resolution_inv_)),
                                        static_cast<int>(std::floor(max(2) * resolution_inv_))}};
            storage->template set<cis::option::tags::array_offset>(static_cast<std::size_t>(min_index[0]),
                                                                   static_cast<std::size_t>(min_index[1]),
                                                                   static_cast<std::size_t>(min_index[2]));
            storage->template set<cis::option::tags::array_size>(static_cast<std::size_t>(max_index[0] - min_index[0] + 1),
                                                                 static_cast<std::size_t>(max_index[1] - min_index[1] + 1),
                                                                 static_cast<std::size_t>(max_index[2] - min_index[2] + 1));
        }

        bundle_storage_->template set<cis::option::tags::array_offset>(static_cast<std::size_t>(min_index_[0]),
                                                                       static_cast<std::size_t>(min_index_[1]),
                                                                       static_cast<std::size_t>(min_index_[2]));
        bundle_storage_->template set<cis::option::tags::array_size>(static_cast<std::size_t>(max_index_[0] - min_index_[0] + 1),
                                                                     static_cast<std::size_t>(max_index_[1] - min_index_[1] + 1),
                                                                     static_cast<std::size_t>(max_index_[2] - min_index_[2] + 1));
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
        return w_T_m_;
    }

    inline void add(const point_t &p)
    {
        lock_t l(bundle_storage_mutex_);
        const index_t bi = toBundleIndex(p);
        distribution_bundle_t *bundle = bundle_storage_->get(bi);
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
    const size_t                                    size_;
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

    inline index_t toIndex(const point_t &p_w,
                           const point_t &off = point_t()) const
    {
        const point_t p_m = (m_T_w_ * p_w) + off;
<<<<<<< HEAD
        return {{static_cast<int>(cslibs_math::common::floor(p_m(0) * resolution_inv_ + 0.5)),
                 static_cast<int>(cslibs_math::common::floor(p_m(1) * resolution_inv_ + 0.5)),
                 static_cast<int>(cslibs_math::common::floor(p_m(2) * resolution_inv_ + 0.5))}};
=======
        return {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * resolution_inv_))}};
>>>>>>> tmp
}

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
<<<<<<< HEAD
        return {{static_cast<int>(cslibs_math::common::floor(p_m(0) * bundle_resolution_inv_ + 0.5)),
                 static_cast<int>(cslibs_math::common::floor(p_m(1) * bundle_resolution_inv_ + 0.5)),
                 static_cast<int>(cslibs_math::common::floor(p_m(2) * bundle_resolution_inv_ + 0.5))}};
=======
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
>>>>>>> tmp

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



#endif // CSLIBS_NDT_3D_STATIC_GRIDMAP_HPP
