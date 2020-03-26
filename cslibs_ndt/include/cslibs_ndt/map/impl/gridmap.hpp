#ifndef CSLIBS_NDT_MAP_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_GRIDMAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>
#include <cslibs_ndt/common/distribution.hpp>

#include <cslibs_ndt/utility/bilinear_interpolation.hpp>

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 Map<option_t,Dim,Distribution,T,backend_t> :
        public GenericMap<option_t,Dim,Distribution,T,backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Map<option_t,Dim,Distribution,T,backend_t>>;

    using ConstPtr = std::shared_ptr<const Map<option_t,Dim,Distribution,T,backend_t>>;
    using Ptr      = std::shared_ptr<Map<option_t,Dim,Distribution,T,backend_t>>;

    using base_t = GenericMap<option_t,Dim,Distribution,T,backend_t>;
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

    using base_t::GenericMap;
    inline Map(const base_t &other) : base_t(other) { }
    inline Map(base_t &&other) : base_t(other) { }

    inline void insert(const typename pointcloud_t::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        return insert(points->begin(), points->end(), points_origin);
    }

    template<typename iterator_t>
    inline void insert(const iterator_t &points_begin,
                       const iterator_t &points_end,
                       const pose_t &points_origin = pose_t())
    {
        std::map<index_t, typename distribution_t::distribution_t> updates;
        for (auto p = points_begin; p != points_end; ++p) {
            const point_t pw = points_origin * *p;
            if (pw.isNormal()) {
                point_t pm;
                index_t bi;
                if (this->toBundleIndex(pw, pm, bi))
                    updates[bi] += pm;
            }
        }

        for (const auto& pair : updates)
            update(pair.first, pair.second);
    }

    inline T sample(const point_t &p) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sample(pm, i);
    }

    inline T sample(const point_t &p,
                    const index_t &bi) const
    {
        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle  = this->bundle_storage_->get(bi);
        return sample(p, bundle);
    }

    inline T sample(const point_t &p,
                    const distribution_bundle_t *bundle) const
    {
        auto evaluate = [this, &p, &bundle]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += this->div_count * bundle->at(i)->data().sample(p);
            return retval;
        };
        return bundle ? evaluate() : T();
    }

    inline T sampleNonNormalized(const point_t &p) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sampleNonNormalized(pm, i);
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const index_t &bi) const
    {
        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle = this->bundle_storage_->get(bi);
        return sampleNonNormalized(p, bundle);
    }

    inline T sampleNonNormalized(const point_t &p,
                                 const distribution_bundle_t *bundle) const
    {
        auto evaluate = [this, &p, &bundle]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += this->div_count * bundle->at(i)->data().sampleNonNormalized(p);
            return retval;
        };
        return bundle ? evaluate() : T();
    }

    inline T sampleNonNormalizedBilinear(const point_t &p) const
    {
        point_t pm;
        const index_t& i = this->toBundleIndex(p, pm);
        return sampleNonNormalizedBilinear(pm, i);
    }

    inline T sampleNonNormalizedBilinear(const point_t &p,
                                         const index_t &bi) const
    {
        if (!this->valid(bi))
            return T();

        distribution_bundle_t *bundle = this->bundle_storage_->get(bi);
        const auto& weights = utility::get_bilinear_interpolation_weights(bi,p,this->bundle_resolution_inv_);
        return sampleNonNormalizedBilinear(p, weights, bundle);
    }

    inline T sampleNonNormalizedBilinear(const point_t &p,
                                         const std::array<T,Dim> &weights,
                                         const distribution_bundle_t *bundle) const
    {
        auto evaluate = [this, &p, &weights, &bundle]() {
            T retval = T();
            for (std::size_t i=0; i<this->bin_count; ++i)
                retval += utility::to_bilinear_interpolation_weight(weights,i) * bundle->at(i)->data().sampleNonNormalized(p);
            return retval;
        };
        return bundle ? evaluate() : T();
    }

protected:
    virtual inline bool expandDistribution(const distribution_t* d) const override
    {
        return d && d->data().valid();
    }

    inline void update(const index_t &bi,
                       const typename distribution_t::distribution_t &d) const
    {
        const distribution_bundle_t *bundle = this->getAllocate(bi);
        for (std::size_t i=0; i<this->bin_count; ++i)
            bundle->at(i)->data() += d;
    }
};
}
}

#endif // CSLIBS_NDT_MAP_GRIDMAP_HPP
