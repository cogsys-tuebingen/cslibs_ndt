#ifndef CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_2d/linear/pose.hpp>
#include <cslibs_math_2d/linear/point.hpp>

#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/linear/pointcloud.hpp>
#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/array/array.hpp>

#include <cslibs_math_2d/algorithms/bresenham.hpp>
#include <cslibs_math_2d/algorithms/simple_iterator.hpp>
#include <cslibs_gridmaps/utility/inverse_model.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_2d {
namespace static_maps {
class EIGEN_ALIGN16 OccupancyGridmap
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using allocator_t = Eigen::aligned_allocator<Gridmap>;

    using ConstPtr                          = std::shared_ptr<const OccupancyGridmap>;
    using Ptr                               = std::shared_ptr<OccupancyGridmap>;
    using pose_t                            = cslibs_math_2d::Pose2d;
    using transform_t                       = cslibs_math_2d::Transform2d;
    using point_t                           = cslibs_math_2d::Point2d;
    using index_t                           = std::array<int, 2>;
    using size_t                            = std::array<std::size_t, 2>;
    using size_m_t                          = std::array<double, 2>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::OccupancyDistribution<2>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::array::Array>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 4>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 4>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 4>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::array::Array>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;
    using simple_iterator_t                 = cslibs_math_2d::algorithms::SimpleIterator;
    using inverse_sensor_model_t            = cslibs_gridmaps::utility::InverseModel;

    inline OccupancyGridmap(const pose_t &origin,
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
        (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
        min_bundle_index[1] + static_cast<int>(size[1] * 2)}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_[0]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                min_bundle_index[1] / 2);
        for(std::size_t i = 1 ; i < 4 ; ++ i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1);
            storage_[i]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                    min_bundle_index[1] / 2);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0],
                min_bundle_index[1]);
    }

    inline OccupancyGridmap(const double &origin_x,
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
        (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
        min_bundle_index[1] + static_cast<int>(size[1] * 2)}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1]);
        storage_[0]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                min_bundle_index[1] / 2);
        for(std::size_t i = 1 ; i < 4 ; ++ i) {
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1);
            storage_[i]->template set<cis::option::tags::array_offset>(min_bundle_index[0] / 2,
                    min_bundle_index[1] / 2);
        }

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2);
        bundle_storage_->template set<cis::option::tags::array_offset>(min_bundle_index[0],
                min_bundle_index[1]);
    }

    inline OccupancyGridmap(const pose_t &origin,
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
        (size[1] + 1) * resolution}},
        min_bundle_index_(min_bundle_index),
        max_bundle_index_{{min_bundle_index[0] + static_cast<int>(size[0] * 2),
        min_bundle_index[1] + static_cast<int>(size[1] * 2)}},
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline OccupancyGridmap(const OccupancyGridmap &other) :
        resolution_(other.resolution_),
        resolution_inv_(other.resolution_inv_),
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

    inline OccupancyGridmap(OccupancyGridmap &&other) :
        resolution_(other.resolution_),
        resolution_inv_(other.resolution_inv_),
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
        cslibs_math_2d::Transform2d origin = w_T_m_;
        origin.translation() += point_t(min_bundle_index_[0] * bundle_resolution_,
                min_bundle_index_[1] * bundle_resolution_);
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

    template <typename line_iterator_t = simple_iterator_t>
    inline void insert(const point_t &start_p,
                       const point_t &end_p)
    {
        const index_t end_index = toBundleIndex(end_p);
        updateOccupied(end_index, end_p);

        line_iterator_t it(m_T_w_ * start_p, m_T_w_ * end_p, bundle_resolution_);
        while (!it.done()) {
            updateFree({{it.x(), it.y()}});
            ++ it;
        }
    }

    template <typename line_iterator_t = simple_iterator_t>
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
                if(toBundleIndex(pm, bi)) {
                    distribution_t *d = storage.get(bi);
                    (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
                }
            }
        }

        const point_t start_p = m_T_w_ * points_origin.translation();
        storage.traverse([this, &start_p](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;
            updateOccupied(bi, d.getDistribution());

            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), bundle_resolution_);
            const std::size_t n = d.numOccupied();
            while (!it.done()) {
                updateFree({{it.x(), it.y()}}, n);
                ++ it;
            }
        });
    }

    template <typename line_iterator_t = simple_iterator_t>
    inline void insertVisible(const pose_t &origin,
                              const typename cslibs_math::linear::Pointcloud<point_t>::ConstPtr &points,
                              const inverse_sensor_model_t::Ptr &ivm,
                              const inverse_sensor_model_t::Ptr &ivm_visibility)
    {
        if (!ivm || !ivm_visibility) {
            std::cout << "[OccupancyGridmap2D]: Cannot evaluate visibility, using model-free update rule instead!" << std::endl;
            return insert(points, origin);
        }

        const index_t start_bi = toBundleIndex(origin.translation());

        auto occupancy = [this, &ivm](const index_t &bi) {
            const distribution_bundle_t *bundle = getDistributionBundle(bi);
            return 0.25 * (bundle->at(0)->getHandle()->getOccupancy(ivm) +
                           bundle->at(1)->getHandle()->getOccupancy(ivm) +
                           bundle->at(2)->getHandle()->getOccupancy(ivm) +
                           bundle->at(3)->getHandle()->getOccupancy(ivm));
        };
        auto current_visibility = [this, &start_bi, &ivm_visibility, &occupancy](const index_t &bi) {
            const double occlusion_prob =
                    std::min(occupancy({{bi[0] + ((bi[0] > start_bi[0]) ? -1 : 1), bi[1]}}),
                             occupancy({{bi[0], bi[1] + ((bi[1] > start_bi[1]) ? -1 : 1)}}));
            return ivm_visibility->getProbFree() * occlusion_prob +
                    ivm_visibility->getProbOccupied() * (1.0 - occlusion_prob);
        };

        distribution_storage_t storage;
        storage.template set<cis::option::tags::array_size>(size_[0] * 2, size_[1] * 2);
        storage.template set<cis::option::tags::array_offset>(min_bundle_index_[0],
                min_bundle_index_[1]);
        for (const auto &p : *points) {
            const point_t pm = origin * p;
            if (pm.isNormal()) {
                index_t bi;
                if(toBundleIndex(pm, bi)) {
                    distribution_t *d = storage.get(bi);
                    (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
                }
            }
        }

        const point_t start_p = m_T_w_ * origin.translation();
        storage.traverse([this, &ivm_visibility, &start_p, &current_visibility](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;

            const point_t end_p = m_T_w_ * point_t(d.getDistribution()->getMean());
            line_iterator_t it(start_p, end_p, bundle_resolution_);

            const std::size_t n = d.numOccupied();
            double visibility = 1.0;
            while (!it.done()) {
                const index_t bit = {{it.x(), it.y()}};
                if ((visibility *= current_visibility(bit)) < ivm_visibility->getProbPrior())
                    return;

                updateFree(bit, n);
                ++ it;
            }

            if ((visibility *= current_visibility(bi)) >= ivm_visibility->getProbPrior())
                updateOccupied(bi, d.getDistribution());
        });
    }

    inline double sample(const point_t &p,
                         const inverse_sensor_model_t::Ptr &ivm) const
    {
        return sample(p, toBundleIndex(p), ivm);
    }

    inline double sample(const point_t &p,
                         const index_t &bi,
                         const inverse_sensor_model_t::Ptr &ivm) const
    {
        if(!valid(bi))
            return 0.0;

        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d->getHandle();
                return handle->getDistribution() ?
                            handle->getDistribution()->sample(p) * handle->getOccupancy(ivm) : 0.0;
            };
            return d ? do_sample() : 0.0;
        };
        auto evaluate = [&p, &bundle, &sample]() {
            return 0.25 * (sample(bundle->at(0)) +
                           sample(bundle->at(1)) +
                           sample(bundle->at(2)) +
                           sample(bundle->at(3)));
        };
        return bundle ? evaluate() : 0.0;
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const inverse_sensor_model_t::Ptr &ivm) const
    {
        return sampleNonNormalized(p, toBundleIndex(p), ivm);
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const index_t &bi,
                                      const inverse_sensor_model_t::Ptr &ivm) const
    {
        if(!valid(bi))
            return 0.0;

        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = bundle_storage_->get(bi);
        }
        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d->getHandle();
                return handle->getDistribution() ?
                            handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : 0.0;
            };
            return d ? do_sample() : 0.0;
        };
        auto evaluate = [&p, &bundle, &sample]() {
            return 0.25 * (sample(bundle->at(0)) +
                           sample(bundle->at(1)) +
                           sample(bundle->at(2)) +
                           sample(bundle->at(3)));
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
        return {{size_[0] * 2, size_[1] * 2}};
    }

    inline distribution_storage_array_t const & getStorages() const
    {
        return storage_;
    }

    template <typename Fn>
    inline void traverse(const Fn& function) const
    {
        lock_t(bundle_storage_mutex_);
        return bundle_storage_->traverse(function);
    }

    inline void getBundleIndices(std::vector<index_t> &indices) const
    {
        lock_t(bundle_storage_mutex_);
        auto add_index = [&indices](const index_t &i, const distribution_bundle_t &d) {
            indices.emplace_back(i);
        };
        bundle_storage_->traverse(add_index);
    }

    inline std::size_t getByteSize() const
    {
        lock_t(bundle_storage_mutex_);
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

        lock_t l(bundle_storage_mutex_);
        const static int dx[] = {-1, 0, 1 -1, 1,-1, 0, 1};
        const static int dy[] = {-1,-1,-1, 0, 0, 1, 1, 1};
        for(const index_t &bi : bis) {
            const distribution_bundle_t *bundle = bundle_storage_->get(bi);
            bool expand = false;
            expand |= bundle->at(0)->getHandle()->getDistribution()->getN() >= 3;
            expand |= bundle->at(1)->getHandle()->getDistribution()->getN() >= 3;
            expand |= bundle->at(2)->getHandle()->getDistribution()->getN() >= 3;
            expand |= bundle->at(3)->getHandle()->getDistribution()->getN() >= 3;

            if(expand) {
                for(std::size_t i = 0 ; i < 8 ; ++i) {
                    const index_t bni = {{dx[i] + bi[0], dy[i] + bi[1]}};
                    if(valid(bni))
                        getAllocate(bni);
                }
            }
        }
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

    mutable distribution_storage_array_t            storage_;
    mutable mutex_t                                 bundle_storage_mutex_;
    mutable distribution_bundle_storage_ptr_t       bundle_storage_;

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
                distribution_bundle_t b;
                const int divx = cslibs_math::common::div(bi[0], 2);
                const int divy = cslibs_math::common::div(bi[1], 2);
                const int modx = cslibs_math::common::mod(bi[0], 2);
                const int mody = cslibs_math::common::mod(bi[1], 2);

                const index_t storage_0_index = {{divx,        divy}};
                const index_t storage_1_index = {{divx + modx, divy}};        /// shifted to the left
                const index_t storage_2_index = {{divx,        divy + mody}}; /// shifted to the bottom
                const index_t storage_3_index = {{divx + modx, divy + mody}}; /// shifted diagonally

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

    inline void updateFree(const index_t &bi) const
    {
        if(!valid(bi))
            return;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateFree();
        bundle->at(1)->getHandle()->updateFree();
        bundle->at(2)->getHandle()->updateFree();
        bundle->at(3)->getHandle()->updateFree();
    }

    inline void updateFree(const index_t &bi,
                           const std::size_t &n) const
    {
        if(!valid(bi))
            return;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateFree(n);
        bundle->at(1)->getHandle()->updateFree(n);
        bundle->at(2)->getHandle()->updateFree(n);
        bundle->at(3)->getHandle()->updateFree(n);
    }

    inline void updateOccupied(const index_t &bi,
                               const point_t &p) const
    {
        if(!valid(bi))
            return;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateOccupied(p);
        bundle->at(1)->getHandle()->updateOccupied(p);
        bundle->at(2)->getHandle()->updateOccupied(p);
        bundle->at(3)->getHandle()->updateOccupied(p);
    }

    inline void updateOccupied(const index_t &bi,
                               const distribution_t::distribution_ptr_t &d) const
    {
        if(!valid(bi))
            return;

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateOccupied(d);
        bundle->at(1)->getHandle()->updateOccupied(d);
        bundle->at(2)->getHandle()->updateOccupied(d);
        bundle->at(3)->getHandle()->updateOccupied(d);
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                        static_cast<int>(std::floor(p_m(1) * resolution_inv_))}};
    }

    inline bool toBundleIndex(const point_t &p_w,
                              index_t &index) const
    {
        const point_t p_m = m_T_w_ * p_w;
        index = {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                  static_cast<int>(std::floor(p_m(1) * resolution_inv_))}};
        return (index[0] >= min_bundle_index_[0] && index[0] <= max_bundle_index_[0] ) &&
                (index[1] >= min_bundle_index_[1] && index[1] <= max_bundle_index_[1] );
    }

    inline bool valid(const index_t &bi) const
    {
        return (bi[0] >= min_bundle_index_[0] && bi[0] <= max_bundle_index_[0] ) &&
                (bi[1] >= min_bundle_index_[1] && bi[1] <= max_bundle_index_[1] );
    }

    inline void fromIndex(const index_t &i,
                          point_t &p_w) const
    {
        p_w = w_T_m_ * point_t(i[0] * resolution_,
                i[1] * resolution_);
    }
};
}
}

#endif // CSLIBS_NDT_2D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
