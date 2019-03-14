#ifndef CSLIBS_NDT_MAP_WEIGHTED_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_WEIGHTED_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>
#include <cslibs_ndt/common/weighted_occupancy_distribution.hpp>

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class EIGEN_ALIGN16 Map<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t> :
        public GenericMap<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Map<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t>>;

    using ConstPtr = std::shared_ptr<const Map<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t>>;
    using Ptr      = std::shared_ptr<Map<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t>>;

    using base_t = GenericMap<option_t,Dim,WeightedOccupancyDistribution,T,backend_t,dynamic_backend_t>;
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
    using typename base_t::dynamic_distribution_storage_t;

    using inverse_sensor_model_t = cslibs_gridmaps::utility::InverseModel<T>;
    using default_iterator_t     = typename map::traits<Dim,T>::default_iterator_t;

    using base_t::GenericMap;
    inline Map(const base_t &other) : base_t(other) { }
    inline Map(base_t &&other) : base_t(other) { }

    template <typename line_iterator_t = default_iterator_t>
    inline void insert(const point_t &start_p,
                       const point_t &end_p)
    {
        const index_t &end_index = this->toBundleIndex(end_p);
        updateOccupied(end_index, end_p);

        line_iterator_t it(this->m_T_w_ * start_p, this->m_T_w_ * end_p, this->bundle_resolution_);
        while (!it.done()) {
            updateFree(it());
            ++ it;
        }
    }

    template <typename line_iterator_t = default_iterator_t>
    inline void insert(const typename cslibs_math::linear::Pointcloud<point_t>::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        dynamic_distribution_storage_t storage;
        for (const auto &p : *points) {
            const point_t pm = points_origin * p;
            if (pm.isNormal()) {
                const index_t &bi = this->toBundleIndex(pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
            }
        }

        const point_t start_p = this->m_T_w_ * points_origin.translation();
        storage.traverse([this, &start_p](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;
            updateOccupied(bi, d.getDistribution());

            line_iterator_t it(start_p, this->m_T_w_ * point_t(d.getDistribution()->getMean()), this->bundle_resolution_);
            const T w = d.weightOccupied();
            while (!it.done()) {
                updateFree(it(), 1, w); // TODO
                ++ it;
            }
        });
    }

    template <typename line_iterator_t = default_iterator_t>
    inline void insertVisible(const pose_t &origin,
                              const typename cslibs_math::linear::Pointcloud<point_t>::ConstPtr &points,
                              const typename inverse_sensor_model_t::Ptr &ivm,
                              const typename inverse_sensor_model_t::Ptr &ivm_visibility)
    {
        if (!ivm || !ivm_visibility) {
            std::cout << "[WeightedOccupancyGridmap2D]: Cannot evaluate visibility, using model-free update rule instead!" << std::endl;
            return insert(points, origin);
        }

        const index_t start_bi = toBundleIndex(origin.translation());
        auto occupancy = [this, &ivm](const index_t &bi) {
            distribution_bundle_t *bundle = this->getAllocate(bi);
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += this->div_count * bundle->at(i)->getOccupancy(ivm);
            return retval;
        };
        auto current_visibility = [this, &start_bi, &ivm_visibility, &occupancy](const index_t &bi) {
            T occlusion_prob = cslibs_math::utility::traits<T>::One;
            auto generate_occlusion_index = [&bi,&start_bi](const std::size_t& counter) {
                index_t retval = bi;
                retval[counter] += ((bi[counter] > start_bi[counter]) ? -1 : 1);
                return retval;
            };

            for (std::size_t i=0; i<Dim; ++i) {
                const index_t test_index = generate_occlusion_index(i);
                if (this->valid(test_index))
                    occlusion_prob = std::min(occlusion_prob, occupancy(test_index));
            }
            return ivm_visibility->getProbFree() * occlusion_prob +
                   ivm_visibility->getProbOccupied() * (cslibs_math::utility::traits<T>::One - occlusion_prob);
        };

        dynamic_distribution_storage_t storage;
        for (const auto &p : *points) {
            const point_t pm = origin * p;
            if (pm.isNormal()) {
                const index_t &bi = this->toBundleIndex(pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
            }
        }

        const point_t start_p = this->m_T_w_ * origin.translation();
        storage.traverse([this, &ivm_visibility, &start_p, &current_visibility](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;

            const point_t end_p = this->m_T_w_ * point_t(d.getDistribution()->getMean());
            line_iterator_t it(start_p, end_p, this->bundle_resolution_);

            const T ww = d.weightOccupied();
            T visibility = cslibs_math::utility::traits<T>::One;
            while (!it.done()) {
                const index_t bit = it();
                if ((visibility *= current_visibility(bit)) < ivm_visibility->getProbPrior())
                    return;

                updateFree(bit, 1, ww);  // TODO!
                ++ it;
            }

            if ((visibility *= current_visibility(bi)) >= ivm_visibility->getProbPrior())
                updateOccupied(bi, d.getDistribution());
        });
    }

    template <typename line_iterator_t = default_iterator_t>
    inline T getRange(const point_t &start_p,
                      const point_t &end_p,
                      const typename inverse_sensor_model_t::Ptr &ivm,
                      const T &occupied_threshold) const
    {
        if (!ivm)
            throw std::runtime_error("[WeightedOccupancyGridmap]: inverse model not set");

        auto to_bundle_index = [this](const point_t &p) {
            index_t retval;
            for (std::size_t i=0; i<Dim; ++i)
                retval[i] = static_cast<int>(std::floor(p(i) * this->bundle_resolution_inv_));
            return retval;
        };
        const index_t start_index = to_bundle_index(start_p);
        const index_t end_index   = to_bundle_index(end_p);
        line_iterator_t it(start_index, end_index);

        auto occupied = [this, &ivm, &occupied_threshold](const index_t &bi) {
            distribution_bundle_t *bundle = this->bundle_storage_->get(bi);
            auto occupancy = [this, &bundle, &ivm]() {
                T retval = T();
                for (std::size_t i=0; i<this->bin_count; ++i)
                    retval += this->div_count * bundle->at(i)->getOccupancy(ivm);
                return retval;
            };
            return bundle && (occupancy() >= occupied_threshold);
        };

        while (!it.done()) {
            if (occupied(it())){
                auto to_point = [this](const index_t& bi) {
                    point_t retval;
                    for (std::size_t i=0; i<Dim; ++i)
                        retval(i) = static_cast<T>(bi[i]) * this->bundle_resolution_;
                    return retval;
                };
                return (start_p - to_point(it())).length();
            }

            ++ it;
        }

        return (start_p - end_p).length();
    }

    inline T sample(const point_t &p,
                    const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        return sample(p, this->toBundleIndex(p), ivm);
    }

    inline T sample(const point_t &p,
                    const index_t &bi,
                    const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[WeightedOccupancyGridmap]: inverse model not set");

        distribution_bundle_t *bundle = this->bundle_storage_->get(bi);

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

    inline T sampleNonNormalized(const point_t &p,
                                 const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        return sampleNonNormalized(p, this->toBundleIndex(p), ivm);
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const index_t &bi,
                                 const typename inverse_sensor_model_t::Ptr &ivm) const
    {
        if (!ivm)
            throw std::runtime_error("[WeightedOccupancyGridmap]: inverse model not set");

        distribution_bundle_t *bundle = this->bundle_storage_->get(bi);

        auto sample = [&p, &ivm] (const distribution_t *d) {
            auto do_sample = [&p, &ivm, &d]() {
                const auto &handle = d;
                return handle->getDistribution() ?
                            handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : T();
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

protected:
    virtual inline bool expandDistribution(const distribution_t* d) const override
    {
        return d && d->getDistribution() && d->getDistribution()->getSampleCount() > 0;
    }

    inline void updateFree(const index_t &bi) const
    {
        distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<Dim; ++i)
            bundle->at(i)->updateFree();
    }

    inline void updateFree(const index_t     &bi,
                           const std::size_t &n,
                           const T           &w) const
    {
        distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<Dim; ++i)
            bundle->at(i)->updateFree(n, w);
    }

    inline void updateOccupied(const index_t &bi,
                               const point_t &p,
                               const T       &w = cslibs_math::utility::traits<T>::One) const
    {
        distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<Dim; ++i)
            bundle->at(i)->updateOccupied(p, w);
    }

    inline void updateOccupied(const index_t &bi,
                               const typename distribution_t::distribution_ptr_t &d) const
    {
        distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<Dim; ++i)
            bundle->at(i)->updateOccupied(d);
    }
};
}
}

#endif // CSLIBS_NDT_MAP_WEIGHTED_OCCUPANCY_GRIDMAP_HPP
