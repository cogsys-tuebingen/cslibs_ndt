#ifndef CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>
#include <cslibs_ndt/common/occupancy_distribution.hpp>
#include <cslibs_math/statistics/mean.hpp>

#include <cslibs_indexed_storage/operations/clustering/grid_neighborhood.hpp>

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          typename T,
          template <typename, typename, typename...> class backend_t,
          template <typename, typename, typename...> class dynamic_backend_t>
class EIGEN_ALIGN16 Map<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t> :
        public GenericMap<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Map<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t>>;

    using ConstPtr = std::shared_ptr<const Map<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t>>;
    using Ptr      = std::shared_ptr<Map<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t>>;

    using base_t = GenericMap<option_t,Dim,OccupancyDistribution,T,backend_t,dynamic_backend_t>;
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

    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using dist_backend_t = cis::backend::kdtree::KDTree<data_interface_t_, index_interface_t_, options_ts_...>;
    template<typename data_interface_t_, typename index_interface_t_, typename... options_ts_>
    using num_backend_t = cis::backend::simple::UnorderedMap<data_interface_t_, index_interface_t_, options_ts_...>;
    struct Num {
        int n;
        inline Num() : n(0) {}
        inline Num(const int& _n) : n(_n) {}
        inline void merge(const Num& num) {n += num.n;}
    };

    using dist_t = typename distribution_t::distribution_t;
    using dist_storage_t = cis::Storage<dist_t, index_t, dist_backend_t>;
    using num_storage_t = cis::Storage<Num, index_t, num_backend_t>;

    using base_t::GenericMap;
    inline Map(const base_t &other) : base_t(other) { }
    inline Map(base_t &&other) : base_t(other) { }
/*
    template <typename line_iterator_t = default_iterator_t>
    inline void insert(const point_t &start_p,
                       const point_t &end_p)
    {
        point_t end_pm;
        const index_t &end_index = this->toBundleIndex(end_p, end_pm);
        updateOccupied(end_index, end_pm);

        line_iterator_t it(this->m_T_w_ * start_p, end_pm, this->bundle_resolution_);
        while (!it.done()) {
            updateFree(it());
            ++ it;
        }
    }*/

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
        dist_storage_t updates;
        for (auto p = points_begin; p != points_end; ++p) {
            if (p->isNormal()) {
                const point_t pw = points_origin * *p;
                if (pw.isNormal()) {
                    point_t pm;
                    index_t bi;
                    if (this->toBundleIndex(pw, pm, bi)) {
                        dist_t *d = updates.get(bi);
                        d ? (*d += pm) : (updates.insert(bi, dist_t(pm)));
                    }
                }
            }
        }

        num_storage_t updates_free;
        const auto& start = this->m_T_w_ * points_origin.translation();
        updates.traverse([this,&start,&updates_free](const index_t& i, const dist_t& d) {
            updateOccupied(i, d);

            line_iterator_t it(start, point_t(d.getMean()), this->bundle_resolution_);//start_bi, pair.first);
            const std::size_t n = d.getN();
            while (!it.done()) {
                const index_t& bi = it();
                Num *num = updates_free.get(bi);
                num ? (num->n += n) : (updates_free.insert(bi, Num(n)));
                ++it;
            }
        }
        );

       /* updates.traverse([&updates_free](const index_t& i, const dist_t& d) {
            Num* num = updates_free.get(i);
            if (num) num->n = 0;
        });
*/
        updates_free.traverse([this,&updates](const index_t& i, const Num& num) {
            //if (num.n > 0)// && !updates.get(i))
                updateFree(i,num.n);
        });
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
/*
        const index_t start_bi = this->toBundleIndex(points_origin.translation());
        auto occupancy = [this, &ivm](const index_t &bi) {
            const distribution_bundle_t *bundle = this->getAllocate(bi);
            T retval = T();
            if (bundle) {
                for (std::size_t i=0; i<this->bin_count; ++i)
                    retval += this->div_count * bundle->at(i)->getOccupancy(ivm);
            }
            return retval;
        };
        auto current_visibility = [this, &start_bi, &ivm_visibility, &occupancy](const index_t &bi) {
            T occlusion_prob = 1.0;
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
                   ivm_visibility->getProbOccupied() * (1.0 - occlusion_prob);
        };

        dynamic_distribution_storage_t storage;
        for (auto p = points_begin; p != points_end; ++p) {
            const point_t pw = points_origin * *p;
            if (pw.isNormal()) {
                point_t pm;
                const index_t &bi = this->toBundleIndex(pw,pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->updateOccupied(pm);
            }
        }

        const point_t start_p = this->m_T_w_ * points_origin.translation();
        storage.traverse([this, &ivm_visibility, &start_p, &current_visibility](const index_t& bi, const distribution_t &d) {
            if (!d.getDistribution())
                return;

            const point_t end_p = point_t(d.getDistribution()->getMean());
            line_iterator_t it(start_p, end_p, this->bundle_resolution_);

            const std::size_t n = d.numOccupied();
            T visibility = 1.0;
            while (!it.done()) {
                const index_t bit = it();
                if ((visibility *= current_visibility(bit)) < ivm_visibility->getProbPrior())
                    return;

                updateFree(bit, n);
                ++ it;
            }

            if ((visibility *= current_visibility(bi)) >= ivm_visibility->getProbPrior())
                updateOccupied(bi, d.getDistribution());
        });*/
    }

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

protected:
    virtual inline bool expandDistribution(const distribution_t* d) const override
    {
        return d && d->getDistribution() && d->getDistribution()->valid();
    }

    inline void updateFree(const index_t     &bi,
                           const std::size_t &n) const
    {
        distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<this->bin_count; ++i)
            bundle->at(i)->updateFree(n);
    }

    inline void updateOccupied(const index_t &bi,
                               const typename distribution_t::distribution_t &d) const
    {
        distribution_bundle_t* bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<this->bin_count; ++i) {
            bundle->at(i)->updateOccupied(d);
        }
    }
};
}
}

#endif // CSLIBS_NDT_MAP_OCCUPANCY_GRIDMAP_HPP
