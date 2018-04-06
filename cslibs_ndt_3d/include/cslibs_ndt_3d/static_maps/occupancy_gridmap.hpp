#ifndef CSLIBS_NDT_3D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_3D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP

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
#include <cslibs_indexed_storage/operations/clustering/grid_neighborhood.hpp>

#include <cslibs_math_3d/algorithms/bresenham.hpp>
#include <cslibs_math_3d/algorithms/simple_iterator.hpp>

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
    using simple_iterator_t                 = cslibs_math_3d::algorithms::SimpleIterator;
    using inverse_sensor_model_t            = cslibs_gridmaps::utility::InverseModel;

    OccupancyGridmap(const pose_t &origin,
                     const double &resolution,
                     const size_t &size) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
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
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1], size[2]);
        for(std::size_t i = 1 ; i < 8 ; ++ i)
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1, size[2] + 1);

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2, size[2] * 2);
    }

    OccupancyGridmap(const double &origin_x,
                     const double &origin_y,
                     const double &origin_phi,
                     const double &resolution,
                     const size_t &size) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
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
        storage_[0]->template set<cis::option::tags::array_size>(size[0], size[1], size[2]);
        for(std::size_t i = 1 ; i < 8 ; ++ i)
            storage_[i]->template set<cis::option::tags::array_size>(size[0] + 1, size[1] + 1, size[2] + 1);

        bundle_storage_->template set<cis::option::tags::array_size>(size[0] * 2, size[1] * 2, size[2] * 2);
    }

    OccupancyGridmap(const pose_t &origin,
                     const double &resolution,
                     const size_t &size,
                     const std::shared_ptr<distribution_bundle_storage_t> &bundles,
                     const distribution_storage_array_t                   &storage) :
        resolution_(resolution),
        resolution_inv_(1.0 / resolution_),
        bundle_resolution_(0.5 * resolution_),
        bundle_resolution_inv_(1.0 / bundle_resolution_),
        w_T_m_(origin),
        m_T_w_(w_T_m_.inverse()),
        size_(size),
        storage_(storage),
        bundle_storage_(bundles)
    {
    }

    inline pose_t getOrigin() const
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
            updateFree({{it.x(), it.y(), it.z()}});
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

            line_iterator_t it(start_p, m_T_w_ * point_t(d.getDistribution()->getMean()), bundle_resolution_);
            const std::size_t n = d.numOccupied();
            while (!it.done()) {
                updateFree({{it.x(), it.y(), it.z()}}, n);
                ++ it;
            }
        });
    }

    template <typename line_iterator_t = simple_iterator_t>
    inline void insertVisible(const pose_t &origin,
                              const typename cslibs_math::linear::Pointcloud<point_t>::Ptr &points,
                              const inverse_sensor_model_t::Ptr &ivm,
                              const inverse_sensor_model_t::Ptr &ivm_visibility)
    {
        if (!ivm || !ivm_visibility) {
            std::cout << "[OccupancyGridmap3D]: Cannot evaluate visibility, using model-free update rule instead!" << std::endl;
            return insert(origin, points);
        }

        const index_t start_bi = toBundleIndex(origin.translation());
        auto occupancy = [this, &ivm](const index_t &bi) {
            const distribution_bundle_t *bundle = getDistributionBundle(bi);
            return 0.125 * (bundle->at(0)->getHandle()->getOccupancy(ivm) +
                            bundle->at(1)->getHandle()->getOccupancy(ivm) +
                            bundle->at(2)->getHandle()->getOccupancy(ivm) +
                            bundle->at(3)->getHandle()->getOccupancy(ivm) +
                            bundle->at(4)->getHandle()->getOccupancy(ivm) +
                            bundle->at(5)->getHandle()->getOccupancy(ivm) +
                            bundle->at(6)->getHandle()->getOccupancy(ivm) +
                            bundle->at(7)->getHandle()->getOccupancy(ivm));
        };
        auto current_visibility = [this, &start_bi, &ivm_visibility, &occupancy](const index_t &bi) {
            const double occlusion_prob =
                    std::min(occupancy({{bi[0] + ((bi[0] > start_bi[0]) ? -1 : 1), bi[1], bi[2]}}),
                             std::min(occupancy({{bi[0], bi[1] + ((bi[1] > start_bi[1]) ? -1 : 1), bi[2]}}),
                                      occupancy({{bi[0], bi[1], bi[2] + ((bi[2] > start_bi[2]) ? -1 : 1)}})));
            return ivm_visibility->getProbFree() * occlusion_prob +
                   ivm_visibility->getProbOccupied() * (1.0 - occlusion_prob);
        };

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
        storage.traverse([this, &ivm_visibility, &start_p, &current_visibility](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;

            const point_t end_p = m_T_w_ * point_t(d.getDistribution()->getMean());
            line_iterator_t it(start_p, end_p, bundle_resolution_);

            const std::size_t n = d.numOccupied();
            double visibility = 1.0;
            while (!it.done()) {
                const index_t bit = {{it.x(), it.y(), it.z()}};
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
        const index_t bi = toBundleIndex(p);
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
        auto evaluate = [&sample, &bundle]() {
            return 0.125 * (sample(bundle->at(0)) +
                            sample(bundle->at(1)) +
                            sample(bundle->at(2)) +
                            sample(bundle->at(3)) +
                            sample(bundle->at(4)) +
                            sample(bundle->at(5)) +
                            sample(bundle->at(6)) +
                            sample(bundle->at(7)));
        };

        return bundle ? evaluate() : 0.0;
    }

    inline double sampleNonNormalized(const point_t &p,
                                      const inverse_sensor_model_t::Ptr &ivm) const
    {
        const index_t bi = toBundleIndex(p);
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
        auto evaluate = [&sample, &bundle]() {
            return 0.125 * (sample(bundle->at(0)) +
                            sample(bundle->at(1)) +
                            sample(bundle->at(2)) +
                            sample(bundle->at(3)) +
                            sample(bundle->at(4)) +
                            sample(bundle->at(5)) +
                            sample(bundle->at(6)) +
                            sample(bundle->at(7)));
        };

        return bundle ? evaluate() : 0.0;
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

protected:
    const double                                    resolution_;
    const double                                    resolution_inv_;
    const double                                    bundle_resolution_;
    const double                                    bundle_resolution_inv_;
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
        lock_t(storage_mutex_);
        distribution_t *d = s->get(i);
        return d ? d : &(s->insert(i, distribution_t()));
    }

    inline distribution_bundle_t *getAllocate(const index_t &bi) const
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
                return &(bundle_storage_->insert(bi, b));
            };
            return bundle ? bundle : allocate_bundle();
        };

        using neighborhood_t = cis::operations::clustering::GridNeighborhoodStatic<std::tuple_size<index_t>::value, 3>;
        static constexpr neighborhood_t grid{};
        grid.visit([this, &get_allocate, &bi](neighborhood_t::offset_t o) {
            const index_t bi_o({{bi[0]+o[0], bi[1]+o[1], bi[2]+o[2]}});
            if (bi_o[0] >= 0 && bi_o[1] >= 0 && bi_o[2] >=0 &&
                    bi_o[0] < (2 * static_cast<int>(size_[0])) &&
                    bi_o[1] < (2 * static_cast<int>(size_[1])) &&
                    bi_o[2] < (2 * static_cast<int>(size_[2])))
                get_allocate(bi_o);
        });

        return get_allocate(bi);
    }

    inline void updateFree(const index_t &bi) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateFree();
        bundle->at(1)->getHandle()->updateFree();
        bundle->at(2)->getHandle()->updateFree();
        bundle->at(3)->getHandle()->updateFree();
        bundle->at(4)->getHandle()->updateFree();
        bundle->at(5)->getHandle()->updateFree();
        bundle->at(6)->getHandle()->updateFree();
        bundle->at(7)->getHandle()->updateFree();
    }

    inline void updateFree(const index_t &bi,
                           const std::size_t &n) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateFree(n);
        bundle->at(1)->getHandle()->updateFree(n);
        bundle->at(2)->getHandle()->updateFree(n);
        bundle->at(3)->getHandle()->updateFree(n);
        bundle->at(4)->getHandle()->updateFree(n);
        bundle->at(5)->getHandle()->updateFree(n);
        bundle->at(6)->getHandle()->updateFree(n);
        bundle->at(7)->getHandle()->updateFree(n);
    }

    inline void updateOccupied(const index_t &bi,
                               const point_t &p) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateOccupied(p);
        bundle->at(1)->getHandle()->updateOccupied(p);
        bundle->at(2)->getHandle()->updateOccupied(p);
        bundle->at(3)->getHandle()->updateOccupied(p);
        bundle->at(4)->getHandle()->updateOccupied(p);
        bundle->at(5)->getHandle()->updateOccupied(p);
        bundle->at(6)->getHandle()->updateOccupied(p);
        bundle->at(7)->getHandle()->updateOccupied(p);
    }

    inline void updateOccupied(const index_t &bi,
                               const distribution_t::distribution_ptr_t &d) const
    {
        distribution_bundle_t *bundle;
        {
            lock_t(bundle_storage_mutex_);
            bundle = getAllocate(bi);
        }
        bundle->at(0)->getHandle()->updateOccupied(d);
        bundle->at(1)->getHandle()->updateOccupied(d);
        bundle->at(2)->getHandle()->updateOccupied(d);
        bundle->at(3)->getHandle()->updateOccupied(d);
        bundle->at(4)->getHandle()->updateOccupied(d);
        bundle->at(5)->getHandle()->updateOccupied(d);
        bundle->at(6)->getHandle()->updateOccupied(d);
        bundle->at(7)->getHandle()->updateOccupied(d);
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

#endif // CSLIBS_NDT_3D_STATIC_MAPS_OCCUPANCY_GRIDMAP_HPP
