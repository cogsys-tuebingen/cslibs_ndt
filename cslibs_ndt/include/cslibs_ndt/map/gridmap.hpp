#ifndef CSLIBS_NDT_MAP_GRIDMAP_HPP
#define CSLIBS_NDT_MAP_GRIDMAP_HPP

#include <cslibs_ndt/map/generic_map.hpp>
#include <cslibs_ndt/common/distribution.hpp>

namespace cslibs_ndt {
namespace map {
template <tags::option option_t,
          std::size_t Dim,
          typename T,
          template <typename, typename, typename...> class backend_t>
class EIGEN_ALIGN16 Gridmap : public GenericMap<option_t,Dim,Distribution,T,backend_t>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using allocator_t = Eigen::aligned_allocator<Gridmap<option_t,Dim,T,backend_t>>;

    using ConstPtr = std::shared_ptr<const Gridmap<option_t,Dim,T,backend_t>>;
    using Ptr      = std::shared_ptr<Gridmap<option_t,Dim,T,backend_t>>;

    inline void insert(const point_t &p)
    {
        index_t bi;
        if (!toBundleIndex(p, bi))
            return;

        distribution_bundle_t *bundle = getAllocate(bi);
        for (std::size_t i=0; i<Dim; ++i)
            bundle->at(i)->data().add(p);
    }

    inline void insert(const typename pointcloud_t::ConstPtr &points,
                       const pose_t &points_origin = pose_t())
    {
        distribution_storage_t storage;
        allocateStorage(storage);

        for (const auto &p : *points) {
            const point_t pm = points_origin * p;
            if (pm.isNormal()) {
                const index_t &bi = toBundleIndex(pm);
                distribution_t *d = storage.get(bi);
                (d ? d : &storage.insert(bi, distribution_t()))->data().add(pm);
            }
        }

        storage.traverse([this](const index_t& bi, const distribution_t &d) {
            distribution_bundle_t *bundle = getAllocate(bi);
            for (std::size_t i=0; i<Dim; ++i)
                bundle->at(i)->data() += d.data();
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
            T retval = 0.0;
            for (std::size_t i=0; i<Dim; ++i)
                retval += 0.25 * bundle->at(i)->data().sample(p);
            return retval;
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
            T retval = 0.0;
            for (std::size_t i=0; i<Dim; ++i)
                retval += 0.25 * bundle->at(i)->data().sampleNonNormalized(p);
            return retval;
        };
        return bundle ? evaluate() : 0.0;
    }
};
}
}

#endif // CSLIBS_NDT_MAP_GRIDMAP_HPP
