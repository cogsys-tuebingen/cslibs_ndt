#ifndef CSLIBS_NDT_3D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP

#include <array>
#include <vector>
#include <cmath>
#include <memory>

#include <cslibs_math_3d/linear/pose.hpp>
#include <cslibs_math_3d/linear/point.hpp>

#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_ndt/common/bundle.hpp>

#include <cslibs_math/linear/pointcloud.hpp>
#include <cslibs_math/common/array.hpp>
#include <cslibs_math/common/div.hpp>
#include <cslibs_math/common/mod.hpp>

#include <cslibs_indexed_storage/storage.hpp>
#include <cslibs_indexed_storage/backend/kdtree/kdtree.hpp>

#include <cslibs_math_3d/algorithms/bresenham.hpp>
#include <cslibs_math_3d/algorithms/simple_iterator.hpp>

#include <unordered_map>
namespace cis = cslibs_indexed_storage;

namespace cslibs_ndt_3d {
namespace dynamic_maps {
class OccupancyGridmap
{
public:
    using Ptr                               = std::shared_ptr<OccupancyGridmap>;
    using pose_t                            = cslibs_math_3d::Pose3d;
    using transform_t                       = cslibs_math_3d::Transform3d;
    using point_t                           = cslibs_math_3d::Point3d;
    using index_t                           = std::array<int, 3>;
    using mutex_t                           = std::mutex;
    using lock_t                            = std::unique_lock<mutex_t>;
    using distribution_t                    = cslibs_ndt::OccupancyDistribution<3>;
    using distribution_storage_t            = cis::Storage<distribution_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_storage_ptr_t        = std::shared_ptr<distribution_storage_t>;
    using distribution_storage_array_t      = std::array<distribution_storage_ptr_t, 8>;
    using distribution_bundle_t             = cslibs_ndt::Bundle<distribution_t*, 8>;
    using distribution_const_bundle_t       = cslibs_ndt::Bundle<const distribution_t*, 8>;
    using distribution_bundle_storage_t     = cis::Storage<distribution_bundle_t, index_t, cis::backend::kdtree::KDTree>;
    using distribution_bundle_storage_ptr_t = std::shared_ptr<distribution_bundle_storage_t>;
    using line_iterator_t                   = cslibs_math_3d::algorithms::SimpleIterator;//Bresenham;

    OccupancyGridmap(const pose_t &origin,
                     const double  resolution) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        bundle_resolution_2_(0.25 * bundle_resolution_ * bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        min_index_{{std::numeric_limits<int>::max(), std::numeric_limits<int>::max()}},
        max_index_{{std::numeric_limits<int>::min(), std::numeric_limits<int>::min()}},
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
        cslibs_math_3d::Transform3d origin = w_T_m_;
        origin.translation() = getMin();
        return origin;
    }

    inline pose_t getInitialOrigin() const
    {
        return w_T_m_;
    }

    inline void add(const point_t &start_p,
                    const point_t &end_p)
    {
        const index_t &end_index = toBundleIndex(end_p);
        line_iterator_t it(m_T_w_ * start_p, m_T_w_ * end_p, bundle_resolution_);

        while (!it.done()) {
            const index_t bi = {{it.x(), it.y(), it.z()}};
            (it.length2() > bundle_resolution_2_) ?
                        updateFree(bi) :
                        updateOccupied(bi, end_p);
            ++ it;
        }
        updateOccupied(end_index, end_p);
    }

    inline void add(const point_t &start_p,
                    const point_t &end_p,
                    index_t       &end_index)
    {
        end_index = toBundleIndex(end_p);
        line_iterator_t it(m_T_w_ * start_p, m_T_w_ * end_p, bundle_resolution_);

        while (!it.done()) {
            const index_t bi = {{it.x(), it.y(), it.z()}};
            (it.length2() > bundle_resolution_2_) ?
                        updateFree(bi) :
                        updateOccupied(bi, end_p);
            ++ it;
        }
        updateOccupied(end_index, end_p);
    }

    inline void insert(const pose_t &origin,
                       const typename cslibs_math::linear::Pointcloud<point_t>::Ptr &points)
    {
        //distribution_storage_t storage;
        std::unordered_map<index_t, int> map;
        for (const auto &p : *points) {
            const point_t pm = origin * p;
            if (pm.isNormal()) {
                const index_t bi = toBundleIndex(pm);
                updateOccupied(bi, pm);
                ++map[bi];

                //const index_t &bi = toBundleIndex(pm);
                //distribution_t *d = storage.get(bi);
                //(d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);

            }
        }

        const index_t start_index = toBundleIndex(origin.translation());
        for (const auto& entry : map) {
            line_iterator_t it(start_index, entry.first);
            while (!it.done()) {
                updateFree({{it.x(), it.y(), it.z()}}, entry.second);
                ++ it;
            }
        }
/*/
        const point_t start_p = m_T_w_ * origin.translation();
        storage.traverse([this, &start_p](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;/*
            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), bundle_resolution_);//start_index, bi, 1.0);

            while (!it.done()) {
                const index_t bj = {{it.x(), it.y(), it.z()}};
                (it.length2() > bundle_resolution_2_) ?
                            updateFree(bj, d.numOccupied()) :
                            updateOccupied(bj, d.getDistribution());
                ++ it;
            }* /

            updateOccupied(bi, d.getDistribution());

            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), bundle_resolution_);//start_index, bi, 1.0);
            const std::size_t n = d.numOccupied();
            while (!it.done()) {
                updateFree({{it.x(), it.y(), it.z()}}, n);
                ++ it;
            }//* /
        });//*/
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
    const double                                    bundle_resolution_2_;
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
            const int divz = cslibs_math::common::div<int>(bi[2], 2);
            const int modx = cslibs_math::common::mod<int>(bi[0], 2);
            const int mody = cslibs_math::common::mod<int>(bi[1], 2);
            const int modz = cslibs_math::common::mod<int>(bi[2], 2);

            const index_t storage_index_0 = {{divx,        divy,        divz}};
            const index_t storage_index_1 = {{divx + modx, divy,        divz}};
            const index_t storage_index_2 = {{divx,        divy + mody, divz}};
            const index_t storage_index_3 = {{divx + modx, divy + mody, divz}};
            const index_t storage_index_4 = {{divx,        divy,        divz + modz}};
            const index_t storage_index_5 = {{divx + modx, divy,        divz + modz}};
            const index_t storage_index_6 = {{divx,        divy + mody, divz + modz}};
            const index_t storage_index_7 = {{divx + modx, divy + mody, divz + modz}};

            b[0] = getAllocate(storage_[0], storage_index_0);
            b[1] = getAllocate(storage_[1], storage_index_1);
            b[2] = getAllocate(storage_[2], storage_index_2);
            b[3] = getAllocate(storage_[3], storage_index_3);
            b[4] = getAllocate(storage_[4], storage_index_4);
            b[5] = getAllocate(storage_[5], storage_index_5);
            b[6] = getAllocate(storage_[6], storage_index_6);
            b[7] = getAllocate(storage_[7], storage_index_7);

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
        bundle->at(4)->updateFree();
        bundle->at(5)->updateFree();
        bundle->at(6)->updateFree();
        bundle->at(7)->updateFree();
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
        bundle->at(4)->updateFree(n);
        bundle->at(5)->updateFree(n);
        bundle->at(6)->updateFree(n);
        bundle->at(7)->updateFree(n);
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
        bundle->at(4)->updateOccupied(d);
        bundle->at(5)->updateOccupied(d);
        bundle->at(6)->updateOccupied(d);
        bundle->at(7)->updateOccupied(d);
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
                 static_cast<int>(std::floor(p_m(1) * bundle_resolution_inv_)),
                 static_cast<int>(std::floor(p_m(2) * bundle_resolution_inv_))}};
    }
};
}
}

#endif // CSLIBS_NDT_3D_DYNAMIC_OCCUPANCY_GRIDMAP_HPP
