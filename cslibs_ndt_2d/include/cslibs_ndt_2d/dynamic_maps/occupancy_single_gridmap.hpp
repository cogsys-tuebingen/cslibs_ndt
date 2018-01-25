#ifndef CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_SINGLE_GRIDMAP_HPP
#define CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_SINGLE_GRIDMAP_HPP

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
#include <cslibs_math_2d/algorithms/efla_iterator.hpp>

namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_2d {
namespace dynamic_maps {
class OccupancySingleGridmap
{
public:
    using Ptr                               = std::shared_ptr<OccupancySingleGridmap>;
    using pose_t                            = cslibs_math_2d::Pose2d;
    using transform_t                       = cslibs_math_2d::Transform2d;
    using point_t                           = cslibs_math_2d::Point2d;
    using index_t                           = std::array<int, 2>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::OccupancyDistribution<2>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using simple_iterator_t                 = cslibs_math_2d::algorithms::SimpleIterator;

    OccupancySingleGridmap(const pose_t &origin,
                           const double  resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        storage_(new distribution_storage_t)
    {
    }

    OccupancySingleGridmap(const pose_t &origin,
                           const double  resolution,
                           const distribution_storage_ptr_t                     &storage,
                           const index_t                                        &min_index,
                           const index_t                                        &max_index) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_(min_index),
        max_index_(max_index),
        storage_(storage)
    {
    }

    OccupancySingleGridmap(const double origin_x,
                           const double origin_y,
                           const double origin_phi,
                           const double resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        w_T_m_(origin_x, origin_y, origin_phi),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
        storage_(new distribution_storage_t)
    {
    }

    inline point_t getMin() const
    {
        lock_t l(storage_mutex_);
        return point_t(min_index_[0] * resolution_,
                min_index_[1] * resolution_);
    }

    inline point_t getMax() const
    {
        lock_t l(storage_mutex_);
        return point_t((max_index_[0] + 1) * resolution_,
                (max_index_[1] + 1) * resolution_);
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
        const index_t &end_index = toIndex(end_p);
        updateOccupied(end_index, end_p);

        line_iterator_t it(m_T_w_ * start_p, m_T_w_ * end_p, resolution_);
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
                const index_t &i = toIndex(pm);
                distribution_t *d = storage.get(i);
                (d ? d : &storage.insert(i, distribution_t()))->updateOccupied(pm);
            }
        }

        const point_t start_p = m_T_w_ * origin.translation();
        storage.traverse([this, &start_p](const index_t& i, const distribution_t &d) {
            if (!d.getDistribution())
                return;
            updateOccupied(i, d.getDistribution());

            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), resolution_);//start_index, i, 1.0);
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
            throw std::runtime_error("[OccupancySingleGridmap]: inverse model not set");

        const index_t start_index = {{static_cast<int>(std::floor(start_p(0) * resolution_inv_)),
                                      static_cast<int>(std::floor(start_p(1) * resolution_inv_))}};
        const index_t end_index   = {{static_cast<int>(std::floor(end_p(0)   * resolution_inv_)),
                                      static_cast<int>(std::floor(end_p(1)   * resolution_inv_))}};
        line_iterator_t it(start_index, end_index);

        auto occupied = [this, &inverse_model, &occupied_threshold](const index_t &i) {
            distribution_t *d;
            {
                lock_t(storage_mutex_);
                d = getAllocate(i);
            }
            return d && d->getOccupancy(inverse_model) >= occupied_threshold;
        };

        while (!it.done()) {
            if (occupied({{it.x(), it.y()}}))
                return (start_p - point_t(it.x() * resolution_, it.y() * resolution_)).length();

            ++ it;
        }

        return (start_p - end_p).length();
    }

    inline double sample(const point_t &p,
                         const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        const index_t i = toIndex(p);
        return sample(p, i, inverse_model);
    }

    inline double sample(const point_t &p,
                         const index_t &i,
                         const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("[OccupancySingleGridmap]: inverse model not set");

        distribution_t *d;
        {
            lock_t(storage_mutex_);
            d = getAllocate(i);
        }
        auto sample = [&p, &inverse_model, &d]() {
            return (d && d->getDistribution()) ? d->getDistribution()->sample(p) *
                                                 d->getOccupancy(inverse_model) : 0.0;
        };
        return sample();
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        const index_t i = toIndex(p);
        return sampleNonNormalized(p, i, inverse_model);
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const index_t &i,
                                      const cslibs_gridmaps::utility::InverseModel::Ptr &inverse_model) const
    {
        if (!inverse_model)
            throw std::runtime_error("[OccupancySingleGridmap]: inverse model not set");

        distribution_t *d;
        {
            lock_t(storage_mutex_);
            d = getAllocate(i);
        }
        auto sampleNonNormalized = [&p, &inverse_model, &d]() {
            return (d && d->getDistribution()) ? d->getDistribution()->sampleNonNormalized(p) *
                                                 d->getOccupancy(inverse_model) : 0.0;
        };
        return sampleNonNormalized();
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

    inline const distribution_t* getDistrbution(const index_t &i) const
    {
        return getAllocate(i);
    }

    inline distribution_t* getDistrbution(const index_t &i)
    {
        return getAllocate(i);
    }

    inline double getResolution() const
    {
        return resolution_;
    }

    inline double getHeight() const
    {
        return (max_index_[1] - min_index_[1] + 1) * resolution_;
    }

    inline double getWidth() const
    {
        return (max_index_[0] - min_index_[0] + 1) * resolution_;
    }

    inline distribution_storage_ptr_t const & getStorage() const
    {
        return storage_;
    }

    inline void getIndices(std::vector<index_t> &indices) const
    {
        lock_t(storage_mutex_);
        auto add_index = [&indices](const index_t &i, const distribution_t &d) {
            indices.emplace_back(i);
        };
        storage_->traverse(add_index);
    }

private:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const transform_t                               w_T_m_;
    const transform_t                               m_T_w_;

    mutable index_t                                 min_index_;
    mutable index_t                                 max_index_;
    mutable mutex_t                                 storage_mutex_;
    mutable distribution_storage_ptr_t              storage_;

    inline distribution_t* getAllocate(const index_t &i) const
    {
        lock_t l(storage_mutex_);
        distribution_t *d = storage_->get(i);
        return d ? d : &(storage_->insert(i, distribution_t()));
    }

    inline void updateFree(const index_t &i) const
    {
        distribution_t *d = nullptr;
        {
            lock_t(storage_mutex_);
            getAllocate(i);
        }
        d->updateFree();
    }

    inline void updateFree(const index_t &i,
                           const std::size_t &n) const
    {
        distribution_t *d = nullptr;
        {
            lock_t(storage_mutex_);
            getAllocate(i);
        }
        d->updateFree(n);
    }

    inline void updateOccupied(const index_t &i,
                               const point_t &p) const
    {
        distribution_t *d = nullptr;
        {
            lock_t(storage_mutex_);
            getAllocate(i);
        }
        d->updateOccupied(p);
    }

    inline void updateOccupied(const index_t &i,
                               const distribution_t::distribution_ptr_t &di) const
    {
        distribution_t *d = nullptr;
        {
            lock_t(storage_mutex_);
            getAllocate(i);
        }
        d->updateOccupied(di);
    }

    inline void updateIndices(const index_t &i) const
    {
        min_index_ = std::min(min_index_, i);
        max_index_ = std::max(max_index_, i);
    }

    inline index_t toIndex(const point_t &p_w) const
    {
        const point_t p_m = m_T_w_ * p_w;
        return {{static_cast<int>(std::floor(p_m(0) * resolution_inv_)),
                        static_cast<int>(std::floor(p_m(1) * resolution_inv_))}};
    }
};
}
}
#endif // CSLIBS_NDT_2D_DYNAMIC_OCCUPANCY_SINGLE_GRIDMAP_HPP
