#ifndef CSLIBS_NDT_2D_STATIC_MAPS_GRIDMAP_HPP
#define CSLIBS_NDT_2D_STATIC_MAPS_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_2d/linear/pose.hpp>
#include <cslibs_math_2d/linear/point.hpp>

#include <cslibs_ndt/common/distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/linear/pointcloud.hpp>
#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>
#include <cslibs_indexed_storage/operations/clustering/grid_neighborhood.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_2d {
namespace static_maps {
template <typename T>
class EIGEN_ALIGN16 Gridmap
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Gridmap<T>>;

    using ConstPtr                          = std::shared_ptr<const Gridmap<T>>;
    using Ptr                               = std::shared_ptr<Gridmap<T>>;
    using pose_t                            = cslibs_math_2d::Pose2d<T>;
    using transform_t                       = cslibs_math_2d::Transform2d<T>;
    using point_t                           = cslibs_math_2d::Point2d<T>;
    using index_t                           = std::array<int, 2>;
    using size_t                            = std::array<std::size_t, 2>;
    using size_m_t                          = std::array<T, 2>;
    using distribution_t                    = cslibs_ndt::Distribution<T,2>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 4>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 4>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 4>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::array::Array>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;

    inline Gridmap(const pose_t &origin,
                   const T      &resolution,
                   const size_t &size,
                   const index_t &min_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2) - 1,
                           min_bundle_index[1] + static_cast<int>(size[1] * 2) - 1}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_[0]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2, min_bundle_index[1] / 2);
        for(std::size_t i = 1 ; i < 4 ; ++ i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1);
            storage_[i]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2, min_bundle_index[1] / 2);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0], min_bundle_index[1]);
    }

    inline Gridmap(const T &origin_x,
                   const T &origin_y,
                   const T &origin_phi,
                   const T &resolution,
                   const size_t &size,
                   const index_t &min_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2) - 1,
                           min_bundle_index[1] + static_cast<int>(size[1] * 2) - 1}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_[0]->template set<cis::option::tags::array_offset>(cslibs_math::common::div<int>(min_bundle_index[0], 2),
                                                                   cslibs_math::common::div<int>(min_bundle_index[1], 2));
        for(std::size_t i = 1 ; i < 4 ; ++ i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1);
            storage_[i]->template set<cis::option::tags::array_offset>(cslibs_math::common::div<int>(min_bundle_index[0], 2),
                                                                       cslibs_math::common::div<int>(min_bundle_index[1], 2));
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0], min_bundle_index[1]);
    }

    inline Gridmap(const pose_t &origin,
                   const T      &resolution,
                   const size_t &size,
                   const std::shared_ptr<distribution_bundle_storage_t> &bundles,
                   const distribution_storage_array_t                   &storage,
                   const index_t &min_bundle_index) :
        resolution_(resolution),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        size_m_{{(size[0] + 1) * resolution,
                 (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2) - 1,
                           min_bundle_index[1] + static_cast<int>(size[1] * 2) - 1}},
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline Gridmap(const Gridmap &other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(other.w_T_m_),
        m_T_w_(other.m_T_w_),
        size_(other.size_),
        size_m_(other.size_m_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_{{distribution_storage_ptr_t(new distribution_storage_t(*other.storage_[0])),
        distribution_storage_ptr_t(new distribution_storage_t(*other.storage_[1])),
        distribution_storage_ptr_t(new distribution_storage_t(*other.storage_[2])),
        distribution_storage_ptr_t(new distribution_storage_t(*other.storage_[3]))}},
        bundle_storage_(new distribution_bundle_storage_t(*other.bundle_storage_))
    {
    }

    inline Gridmap(Gridmap &&other) :
        resolution_(other.resolution_),
        bundle_resolution_(other.bundle_resolution_),
        bundle_resolution_inv_(other.bundle_resolution_inv_),
        w_T_m_(std::move(other.w_T_m_)),
        m_T_w_(std::move(other.m_T_w_)),
        size_(other.size_),
        size_m_(other.size_m_),
        min_bundle_index_(other.min_bundle_index_),
        max_bundle_index_(other.max_bundle_index_),
        storage_(other.storage_),
        bundle_storage_(other.bundle_storage_)
    {
    }

    /**
     * @brief Get minimum in map coordinates.
     * @return the minimum
     */
    inline point_t getMin() const
    {
        return point_t(min_bundle_index_[0] * bundle_resolution_,
                       min_bundle_index_[1] * bundle_resolution_);
    }

    /**
     * @brief Get maximum in map coordinates.
     * @return the maximum
     */
    inline point_t getMax() const
    {
        return point_t((max_bundle_index_[0] + 1) * bundle_resolution_,
                       (max_bundle_index_[1] + 1) * bundle_resolution_);
    }

    /**
     * @brief Get the origin.
     * @return the origin
     */
    inline pose_t getOrigin() const
    {
        pose_t origin = w_T_m_;
        origin.translation() += getMin();
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
        index_t bi;
        if(!toBundleIndex(p, bi))
            return;

        distribution_bundle_t *bundle = getAllocate(bi);
        bundle->at(0)->data().add(p);
        bundle->at(1)->data().add(p);
        bundle->at(2)->data().add(p);
        bundle->at(3)->data().add(p);
    }

    inline void insert(const typename cslibs_math::linear::Pointcloud<point_t>::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        distribution_storage_t storage;
        storage.template set<cis::option::tags::array_size>(size_[0] * 2, size_[1] * 2);
        storage.template set<cis::option::tags::array_offset>(min_bundle_index_[0],
                                                              min_bundle_index_[1]);
        for (const auto &p : *points) {
            const point_t pm = points_origin * p;
            if (pm.isNormal()) {
                index_t bi;
                if(toBundleIndex(p, bi)) {
                    distribution_t *d = storage.get(bi);
                    (d ? d : &storage.insert(bi, distribution_t()))->data().add(pm);
                }
            }
        }

        storage.traverse([this](const index_t& bi, const distribution_t &d) {
            distribution_bundle_t *bundle = getAllocate(bi);
            bundle->at(0)->data() += d.data();
            bundle->at(1)->data() += d.data();
            bundle->at(2)->data() += d.data();
            bundle->at(3)->data() += d.data();
        });
    }

    inline T sample(const point_t &p) const
    {
        return sample(p, toBundleIndex(p));
    }

    inline T sample(const point_t &p,
                    const index_t &bi) const
    {
        if(!valid(bi))
            return 0.0;

        distribution_bundle_t *bundle  = bundle_storage_->get(bi);
        auto evaluate = [&p, &bundle]() {
            return 0.25 * (bundle->at(0)->data().sample(p) +
                           bundle->at(1)->data().sample(p) +
                           bundle->at(2)->data().sample(p) +
                           bundle->at(3)->data().sample(p));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline T sampleNonNormalized(const point_t &p) const
    {
        return sampleNonNormalized(p, toBundleIndex(p));
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const index_t &bi) const
    {
        if(!valid(bi))
            return 0.0;

        distribution_bundle_t *bundle = bundle_storage_->get(bi);

        auto evaluate = [&p, &bundle]() {
            return 0.25 * (bundle->at(0)->data().sampleNonNormalized(p) +
                           bundle->at(1)->data().sampleNonNormalized(p) +
                           bundle->at(2)->data().sampleNonNormalized(p) +
                           bundle->at(3)->data().sampleNonNormalized(p));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline const distribution_bundle_t* getDistributionBundle(const index_t &bi) const
    {
        return valid(bi) ? getAllocate(bi) : nullptr;
    }

    inline distribution_bundle_t* getDistributionBundle(const index_t &bi)
    {
        return valid(bi) ? getAllocate(bi) : nullptr;
    }

    inline T getBundleResolution() const
    {
        return bundle_resolution_;
    }

    inline T getResolution() const
    {
        return resolution_;
    }

    inline T getHeight() const
    {
        return size_m_[1];
    }

    inline T getWidth() const
    {
        return size_m_[0];
    }

    inline size_t getSize() const
    {
        return size_;
    }

    inline size_t getBundleSize() const
    {
        return {{size_[0] * 2, size_[1] * 2}};
    }

    inline distribution_storage_array_t const & getStorages() const
    {
        return storage_;
    }

    inline index_t getMinBundleIndex() const
    {
        return min_bundle_index_;
    }

    inline index_t getMaxBundleIndex() const
    {
        return max_bundle_index_;
    }

    template <typename Fn>
    inline void traverse(const Fn& function) const
    {
        return bundle_storage_->traverse(function);
    }

    inline void getBundleIndices(std::vector<index_t> &indices) const
    {
        auto add_index = [&indices](const index_t &i, const distribution_bundle_t &d) {
            indices.emplace_back(i);
        };
        bundle_storage_->traverse(add_index);
    }

    inline std::size_t getByteSize() const
    {
        return sizeof(*this) +
                bundle_storage_->byte_size() +
                storage_[0]->byte_size() +
                storage_[1]->byte_size() +
                storage_[2]->byte_size() +
                storage_[3]->byte_size();
    }

    inline virtual bool validate(const pose_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w.translation();
        const point_t min = getMin();
        const point_t max = getMax();
        return p_m(0) >= min(0) && p_m(0) < min(0) &&
               p_m(1) >= min(1) && p_m(1) < max(1);
    }

    inline void allocatePartiallyAllocatedBundles()
    {
        std::vector<index_t> bis;
        getBundleIndices(bis);

        using neighborhood_t = cis::operations::clustering::GridNeighborhoodStatic<std::tuple_size<index_t>::value, 3>;
        static constexpr neighborhood_t grid{};

        for(const index_t &bi : bis) {
            const distribution_bundle_t *bundle  = bundle_storage_->get(bi);
            bool expand =
                    (bundle->at(0)->data().getN() >= 3) ||
                    (bundle->at(1)->data().getN() >= 3) ||
                    (bundle->at(2)->data().getN() >= 3) ||
                    (bundle->at(3)->data().getN() >= 3);

            if(expand) {
                grid.visit([this, &bi](neighborhood_t::offset_t o) {
                    const index_t bni = {{bi[0]+o[0], bi[1]+o[1]}};
                    if (valid(bni))
                        getAllocate(bni);
                });
            }
        }
    }

protected:
    const T                                    resolution_;
    const T                                    bundle_resolution_;
    const T                                    bundle_resolution_inv_;
    const transform_t                          w_T_m_;
    const transform_t                          m_T_w_;
    const size_t                               size_;
    const size_m_t                             size_m_;
    const index_t                              min_bundle_index_;
    const index_t                              max_bundle_index_;

    mutable distribution_storage_array_t       storage_;
    mutable distribution_bundle_storage_ptr_t  bundle_storage_;

    inline distribution_t* getAllocate(const distribution_storage_ptr_t &s,
                                       const index_t &i) const
    {
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
    {
        auto get_allocate = [this](const index_t &bi) {
            distribution_bundle_t *bundle = bundle_storage_->get(bi);

            auto allocate_bundle = [this, &bi]() {
                const int divx = cslibs_math::common::div(bi[0], 2);
                const int divy = cslibs_math::common::div(bi[1], 2);
                const int modx = cslibs_math::common::mod(bi[0], 2);
                const int mody = cslibs_math::common::mod(bi[1], 2);

                const index_t storage_0_index = {{divx,        divy}};
                const index_t storage_1_index = {{divx + modx, divy}};        /// shifted to the left
                const index_t storage_2_index = {{divx,        divy + mody}}; /// shifted to the bottom
                const index_t storage_3_index = {{divx + modx, divy + mody}}; /// shifted diagonally

                distribution_bundle_t b;
                b[0] = getAllocate(storage_[0], storage_0_index);
                b[1] = getAllocate(storage_[1], storage_1_index);
                b[2] = getAllocate(storage_[2], storage_2_index);
                b[3] = getAllocate(storage_[3], storage_3_index);

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
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_))}};
    }

    inline bool toBundleIndex(const point_t &p_w,
                              index_t &index) const
    {
        const point_t p_m = m_T_w_ * p_w;
        index = {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                  static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_))}};
        return (index[0] >= min_bundle_index_[0] && index[0] <= max_bundle_index_[0]) &&
               (index[1] >= min_bundle_index_[1] && index[1] <= max_bundle_index_[1]);
    }

    inline bool valid(const index_t &bi) const
    {
        return (bi[0] >= min_bundle_index_[0] && bi[0] <= max_bundle_index_[0]) &&
               (bi[1] >= min_bundle_index_[1] && bi[1] <= max_bundle_index_[1]);
    }
};
}
}

#endif // CSLIBS_NDT_2D_STATIC_MAPS_GRIDMAP_HPP
