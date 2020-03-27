#ifndef CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>
#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_math/statistics/mean.hpp>

#include <cslibs_indexed_storage/operations/clustering/grid_neighborhood.hpp>
#include <set>

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 Map<option_t,Dim,OccupancyDistribution,T,backend_t> :
        public GenericMap<option_t,Dim,OccupancyDistribution,T,backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Map<option_t,Dim,OccupancyDistribution,T,backend_t>>;

    using ConstPtr = std::shared_ptr<const Map<option_t,Dim,OccupancyDistribution,T,backend_t>>;
    using Ptr      = std::shared_ptr<Map<option_t,Dim,OccupancyDistribution,T,backend_t>>;

    using base_t = GenericMap<option_t,Dim,OccupancyDistribution,T,backend_t>;
    using typename base_t::pose_t;
    using typename base_t::transform_t;
    using typename base_t::point_t;
    using typename base_t::pointcloud_t;
    using typename base_t::index_t;
    using typename base_t::index_list_t;
    using typename base_t::distribution_t;
    using typename base_t::distribution_storage_t;
    using typename base_t::distribution_storage_ptr_t;
    using typename base_t::distribution_storage_array_t;
    using typename base_t::distribution_bundle_t;
    using typename base_t::distribution_const_bundle_t;
    using typename base_t::distribution_bundle_storage_t;
    using typename base_t::distribution_bundle_storage_ptr_t;

    using inverse_sensor_model_t = cslibs_gridmaps::utility::InverseModel<T>;
    using default_iterator_t     = typename map::traits<Dim,T>::default_iterator_t;

    using base_t::GenericMap;
    inline Map(const base_t &other) : base_t(other) { }
    inline Map(base_t &&other) : base_t(other) { }

    template <typename line_iterator_t = default_iterator_t>
    inline void insert(const typename pointcloud_t::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        return insert<line_iterator_t>(points->begin(), points->end(), points_origin);
    }

    template <typename line_iterator_t = default_iterator_t, typename iterator_t>
    inline void insert(const iterator_t &points_begin,
                       const iterator_t &points_end,
                       const pose_t &points_origin = pose_t())
    {
        using dist_t = typename distribution_t::distribution_t;
        std::map<index_t, dist_t> updates;
        for (auto p = points_begin; p != points_end; ++p) {
            if (p->isNormal()) {
                const point_t pw = points_origin * *p;
                if (pw.isNormal()) {
                    point_t pm;
                    index_t bi;
                    if (this->toBundleIndex(pw, pm, bi))
                        updates[bi] += pm;
                }
            }
        }

        std::unordered_map<index_t,std::size_t> updates_free;
        const auto& start = this->m_T_w_ * points_origin.translation();
        for (const auto& pair : updates) {
            const index_t& i = pair.first;
            const dist_t&  d = pair.second;

            const auto& n = d.getN();
            if (this->valid(i))
                updateOccupied(i, d);

            line_iterator_t it(start, point_t(d.getMean()), this->bundle_resolution_);
            while (!it.done()) {
                const index_t& bi = it();
                if (this->valid(bi))
                    updates_free[bi] += n;
                ++it;
            }
        }

        for (const auto& pair : updates)
            updates_free.erase(pair.first);

        for (const auto& pair : updates_free)
            updateFree(pair.first, pair.second);
    }

    template <typename line_iterator_t = default_iterator_t>
    inline void insertVisible(const typename pointcloud_t::ConstPtr &points,
                              const pose_t &points_origin,
                              const typename inverse_sensor_model_t::Ptr &ivm,
                              const typename inverse_sensor_model_t::Ptr &ivm_visibility)
    {
        return insertVisible<line_iterator_t>(points->begin(), points->end(), points_origin, ivm, ivm_visibility);
    }

    template <typename line_iterator_t = default_iterator_t, typename iterator_t>
    inline void insertVisible(const iterator_t &points_begin,
                              const iterator_t &points_end,
                              const pose_t &points_origin,
                              const typename inverse_sensor_model_t::Ptr &ivm,
                              const typename inverse_sensor_model_t::Ptr &ivm_visibility)
    {
        if (!ivm || !ivm_visibility) {
            std::cout << "[OccupancyGridmap]: Cannot evaluate visibility, using model-free update rule instead!" << std::endl;
            return insert<line_iterator_t>(points_begin, points_end, points_origin);
        }

        using dist_t = typename distribution_t::distribution_t;
        std::map<index_t, dist_t> updates;
        for (auto p = points_begin; p != points_end; ++p) {
            if (p->isNormal()) {
                const point_t pw = points_origin * *p;
                if (pw.isNormal()) {
                    point_t pm;
                    index_t bi;
                    if (this->toBundleIndex(pw, pm, bi))
                        updates[bi] += pm;
                }
            }
        }

        std::unordered_map<index_t,std::size_t> updates_free;
        const auto& start = this->m_T_w_ * points_origin.translation();

        const index_t start_bi = this->toBundleIndex(points_origin.translation());
        auto current_visibility = [this, &ivm, &start_bi, &ivm_visibility, &points_origin](const index_t &bi, const point_t& end) {
            auto generate_occlusion_index = [&bi,&start_bi](const std::size_t& counter) {
                index_t retval = bi;
                retval[counter] += ((bi[counter] > start_bi[counter]) ? -1 : 1);
                return retval;
            };
            auto occupancy = [this, &ivm, &points_origin, &end](const index_t &bi) {
                const distribution_bundle_t *bundle = this->get(bi);
                T retval = T(0.);
                if (bundle) {
                    for (std::size_t i=0; i<this->bin_count; ++i) {
                        retval += this->div_count * bundle->at(i)->getOccupancy(ivm);
                    }
                }
                return retval;
            };

            T occlusion_prob = T(1.);
            for (std::size_t i=0; i<Dim; ++i) {
                const index_t test_index = generate_occlusion_index(i);
                if (this->valid(test_index))
                    occlusion_prob = std::min(occlusion_prob, occupancy(test_index));
            }
            return ivm_visibility->getProbFree() * occlusion_prob +
                   ivm_visibility->getProbOccupied() * (T(1.) - occlusion_prob);
        };

        for (const auto& pair : updates) {
            const index_t& i = pair.first;
            const dist_t&  d = pair.second;

            T visibility = T(1.);
            const auto& end = point_t(d.getMean());
            const auto& n = d.getN();

            line_iterator_t it(start, end, this->bundle_resolution_);
            while (!it.done()) {
                const index_t& bi = it();
                if ((visibility *= current_visibility(bi,end)) < ivm_visibility->getProbPrior())
                    return;

                if (this->valid(bi))
                    updates_free[bi] += n;
                ++it;
            }

            if ((visibility *= current_visibility(i,end)) >= ivm_visibility->getProbPrior()) {
                updateOccupied(i, d);
            }
        }

        for (const auto& pair : updates)
            updates_free.erase(pair.first);

        for (const auto& pair : updates_free)
            updateFree(pair.first, pair.second);
    }
/*
    inline T sample(const point_t &p,
                    const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sample(pm, i, ivm);
    }

    inline T sample(const point_t &p,
                    const index_t &bi,
                    const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle = this->bundle_storage_->get(bi);
        return sample(p, bundle, ivm);
    }

    inline T sample(const point_t &p,
                    const distribution_bundle_t* bundle,
                    const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d;
                return handle->getDistribution() ?
                            handle->getDistribution()->sample(p) * handle->getOccupancy(ivm) : T();
            };
            return d ? do_sample() : T();
        };

        auto evaluate = [this, &p, &bundle, &sample]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += this->div_count * sample(bundle->at(i));
            return retval;
        };
        return bundle ? evaluate() : T();
    }
*/
    inline T sampleNonNormalized(const point_t &p,
                                 const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sampleNonNormalized(pm, i, ivm);
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const index_t &bi,
                                 const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle  = this->bundle_storage_->get(bi);
        return sampleNonNormalized(p, bundle, ivm);
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const distribution_bundle_t* bundle,
                                 const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d;
                return handle->getDistribution() ?
                            handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : T(0.0);
            };
            return d ? do_sample() : T();
        };

        auto evaluate = [this, &p, &bundle, &sample]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += this->div_count * sample(bundle->at(i));
            return retval;
        };
        return bundle ? evaluate() : T();
    }

    inline T sampleNonNormalizedBilinear(const point_t &p,
                                         const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sampleNonNormalizedBilinear(pm, i, ivm);
    }

    inline T sampleNonNormalizedBilinear(const point_t &p,
                                         const index_t &bi,
                                         const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle  = this->bundle_storage_->get(bi);
        const auto& weights = utility::get_bilinear_interpolation_weights(bi,p,this->bundle_resolution_inv_);
        return sampleNonNormalizedBilinear(p, weights, bundle, ivm);
    }

    inline T sampleNonNormalizedBilinear(const point_t &p,
                                         const std::array<T,Dim> &weights,
                                         const distribution_bundle_t* bundle,
                                         const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[OccupancyGridMap]: inverse model not set");

        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d;
                return handle->getDistribution() ?
                            handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : T(0.0);
            };
            return d ? do_sample() : T();
        };

        auto evaluate = [this, &p, &weights, &bundle, &sample]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += utility::to_bilinear_interpolation_weight(weights,i) * sample(bundle->at(i));
            return retval;
        };
        return bundle ? evaluate() : T();
    }

protected:
    virtual inline bool expandDistribution(const distribution_t* d) const override
    {
        return d && d->getDistribution() && d->getDistribution()->valid();
    }

    inline void updateFree(const index_t &bi,
                           const std::size_t &n) const
    {
        const distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<this->bin_count; ++i)
            bundle->at(i)->updateFree(n);
    }

    inline void updateOccupied(const index_t &bi,
                               const typename distribution_t::distribution_t &d) const
    {
        const distribution_bundle_t* bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<this->bin_count; ++i)
            bundle->at(i)->updateOccupied(d);
    }
};
}
}

#endif // CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP
