#ifndef CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP

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
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_gridmaps/utility/inverse_model.hpp>

#include <cslibs_math_2d/algorithms/bresenham.hpp>
#include <cslibs_math_2d/algorithms/simple_iterator.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_2d {
namespace dynamic_maps {
class OccupancyGridmap
{
public:
    using Ptr                               = std::shared_ptr<OccupancyGridmap>;
    using pose_t                            = cslibs_math_2d::Pose2d;
    using transform_t                       = cslibs_math_2d::Transform2d;
    using point_t                           = cslibs_math_2d::Point2d;
    using index_t                           = std::array<int, 2>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::OccupancyDistribution<2>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 4>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 4>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 4>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;
    using simple_iterator_t                 = cslibs_math_2d::algorithms::SimpleIterator;

    OccupancyGridmap(const pose_t &origin,
                     const double  resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t),
                 distribution_storage_ptr_t(new distribution_storage_t)}},
        bundle_storage_(new distribution_bundle_storage_t)
    {
    }

    OccupancyGridmap(const double origin_x,
                     const double origin_y,
                     const double origin_phi,
                     const double resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        storage_{{distribution_storage_ptr_t(new distribution_storage_t),
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
                       min_index_[1] * bundle_resolution_);
    }

    inline point_t getMax() const
    {
        lock_t l(bundle_storage_mutex_);
        return point_t((max_index_[0] + 1) * bundle_resolution_,
                       (max_index_[1] + 1) * bundle_resolution_);
    }

    inline pose_t getOrigin() const
    {
        cslibs_math_2d::Transform2d origin = w_T_m_;
        origin.translation() = getMin();
        return origin;
    }

    inline pose_t getInitialOrigin() const
    {
        return w_T_m_;
    }

    template <typename line_iterator_t = simple_iterator_t>
    inline void add(const point_t &start_p,
                    const point_t &end_p)
    {
        const index_t &end_index = toBundleIndex(end_p);
        updateOccupied(end_index, end_p);

        line_iterator_t it(m_T_w_ * start_p, m_T_w_ * end_p, bundle_resolution_);
        while (!it.done()) {
            updateFree({{it.x(), it.y()}});
            ++ it;
        }
    }

    template <typename line_iterator_t = simple_iterator_t>
    inline void insert(const pose_t &origin,
                       const typename cslibs_math::linear::Pointcloud<point_t>::Ptr &points)
    {
        distribution_storage_t storage;
        for (const auto &p : *points) {
            const point_t pm = origin * p;
            if (pm.isNormal()) {
                const index_t &bi = toBundleIndex(pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
            }
        }

        const point_t start_p = m_T_w_ * origin.translation();
        storage.traverse([this, &start_p](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;
            updateOccupied(bi, d.getDistribution());

            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), bundle_resolution_);//start_index, bi, 1.0);
            const std::size_t n = d.numOccupied();
            while (!it.done()) {
                updateFree({{it.x(), it.y()}}, n);
                ++ it;
            }
        });
    }

    template <typename line_iterator_t = simple_iterator_t>
    inline double getRange(const point_t &start_p,
                           const point_t &end_p,
                           const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model,
                           const double &occupied_threshold) const
    {
        if (!inverse_model)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        const index_t start_index = {{static_cast<int>(std::floor(start_p(0) * bundle_resolution_inv_)),
                                      static_cast<int>(std::floor(start_p(1) * bundle_resolution_inv_))}};
        const index_t end_index   = {{static_cast<int>(std::floor(end_p(0) * bundle_resolution_inv_)),
                                      static_cast<int>(std::floor(end_p(1) * bundle_resolution_inv_))}};
        line_iterator_t it(start_index, end_index);

        auto occupied = [this, &inverse_model, &occupied_threshold](const index_t &bi) {
            distribution_bundle_t *bundle;
            {
                lock_t(bundle_storage_mutex_);
                bundle = getAllocate(bi);
            }
            return bundle && (0.25 * ((bundle->at(0)->getOccupancy(inverse_model)) +
                                      (bundle->at(1)->getOccupancy(inverse_model)) +
                                      (bundle->at(2)->getOccupancy(inverse_model)) +
                                      (bundle->at(3)->getOccupancy(inverse_model))) >= occupied_threshold);
        };

        while (!it.done()) {
            if (occupied({{it.x(), it.y()}}))
                return (start_p - point_t(it.x() * bundle_resolution_, it.y() * bundle_resolution_)).length();

            ++ it;
        }

        return (start_p - end_p).length();
    }

    inline double sample(const point_t &p,
                         const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        const index_t bi = toBundleIndex(p);
        return sample(p, bi, inverse_model);
    }

    inline double sample(const point_t &p,
                         const index_t &bi,
                         const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        auto sample = [&p, &inverse_model](const distribution_t *d) {
            return (d && d->getDistribution()) ? d->getDistribution()->sample(p) *
                                                 d->getOccupancy(inverse_model) : 0.0;
        };
        auto evaluate = [&p, &bundle, &sample]() {
            return 0.25 * (sample(bundle->at(0)) +
                           sample(bundle->at(1)) +
                           sample(bundle->at(2)) +
                           sample(bundle->at(3)));
        };
        return evaluate();
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        const index_t bi = toBundleIndex(p);
        return sampleNonNormalized(p, bi, inverse_model);
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const index_t &bi,
                                      const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        auto sampleNonNormalized = [&p, &inverse_model](const distribution_t *d) {
            return (d && d->getDistribution()) ? d->getDistribution()->sampleNonNormalized(p) *
                                                 d->getOccupancy(inverse_model) : 0.0;
        };
        auto evaluate = [&p, &bundle, &sampleNonNormalized]() {
          return 0.25 * (sampleNonNormalized(bundle->at(0)) +
                         sampleNonNormalized(bundle->at(1)) +
                         sampleNonNormalized(bundle->at(2)) +
                         sampleNonNormalized(bundle->at(3)));
        };
        return evaluate();
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
        return (max_index_[1] - min_index_[1] + 1) * bundle_resolution_;
    }

    inline double getWidth() const
    {
        return (max_index_[0] - min_index_[0] + 1) * bundle_resolution_;
    }

    inline distribution_storage_array_t const & getStorages() const
    {
        return storage_;
    }

    inline void getBundleIndices(std::vector<index_t> &indices) const
    {
        lock_t(bundle_storage_mutex_);
        auto add_index = [&indices](const index_t &i, const distribution_bundle_t &d) {
            indices.emplace_back(i);
        };
        bundle_storage_->traverse(add_index);
    }

private:
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
            const int modx = cslibs_math::common::mod<int>(bi[0], 2);
            const int mody = cslibs_math::common::mod<int>(bi[1], 2);

            const index_t storage_index_0 = {{divx,        divy}};
            const index_t storage_index_1 = {{divx + modx, divy}};
            const index_t storage_index_2 = {{divx,        divy + mody}};
            const index_t storage_index_3 = {{divx + modx, divy + mody}};

            b[0] = getAllocate(storage_[0], storage_index_0);
            b[1] = getAllocate(storage_[1], storage_index_1);
            b[2] = getAllocate(storage_[2], storage_index_2);
            b[3] = getAllocate(storage_[3], storage_index_3);

            lock_t(bundle_storage_mutex_);
            updateIndices(bi);
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
    }

    inline void updateFree(const index_t &bi,
                           const std::size_t &n) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->updateFree(n);
        bundle->at(1)->updateFree(n);
        bundle->at(2)->updateFree(n);
        bundle->at(3)->updateFree(n);
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
    }

    inline void updateOccupied(const index_t &bi,
                               const distribution_t::distribution_ptr_t &d) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->updateOccupied(d);
        bundle->at(1)->updateOccupied(d);
        bundle->at(2)->updateOccupied(d);
        bundle->at(3)->updateOccupied(d);
    }

    inline void updateIndices(const index_t &bi) const
    {
        min_index_ = std::min(min_index_, bi);
        max_index_ = std::max(max_index_, bi);
    }

    inline index_t toBundleIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return {{static_cast<int>(std::floor(p_m(0) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_))}};
    }
};
}
}

#endif // CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP
