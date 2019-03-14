#ifndef CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
#define CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_ndt_3d/DistributionArray.h>

namespace cslibs_ndt_3d {
namespace conversion {
template <typename T>
inline Distribution from(const cslibs_math::statistics::Distribution<T, 3, 3> &d,
                         const int &id,
                         const T &prob)
{
    Distribution distr;
    distr.id.data = id;
    for (int i = 0; i < 3; ++ i) {
        distr.mean[i].data          = d.getMean()(i);
        distr.eigen_values[i].data  = d.getEigenValues()(i);
    }
    for (int i = 0; i < 9; ++ i) {
        distr.eigen_vectors[i].data = d.getEigenVectors()(i);
        distr.covariance[i].data    = d.getCovariance()(i);
    }
    distr.prob.data = prob;
    return distr;
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::Gridmap<T>::Ptr &src,
        cslibs_ndt_3d::DistributionArray::Ptr &dst)
{
    if (!src)
        return;
    src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_3d::dynamic_maps::Gridmap<T>;
    using dst_map_t = cslibs_ndt_3d::DistributionArray;
    dst.reset(new dst_map_t());

    using point_t               = typename src_map_t::point_t;
    using distribution_t        = typename src_map_t::distribution_t;
    using distribution_bundle_t = typename src_map_t::distribution_bundle_t;
    auto sample = [](const distribution_t *d,
                     const point_t &p) -> T {
        return d ? d->data().sampleNonNormalized(p) : 0.0;
    };
    auto sample_bundle = [&sample](const distribution_bundle_t &b,
                                   const point_t &p) -> T {
        return 0.125 * (sample(b.at(0), p) +
                        sample(b.at(1), p) +
                        sample(b.at(2), p) +
                        sample(b.at(3), p) +
                        sample(b.at(4), p) +
                        sample(b.at(5), p) +
                        sample(b.at(6), p) +
                        sample(b.at(7), p));
    };

    using index_t = std::array<int, 3>;
    auto process_bundle = [&dst, &sample_bundle](const index_t &bi, const distribution_bundle_t &b) {
        typename distribution_t::distribution_t d;
        for (std::size_t i = 0; i < 8; ++ i)
            d += b.at(i)->data();
        if (d.getN() == 0)
            return;

        dst->data.emplace_back(from(d, b.id(), sample_bundle(b, point_t(d.getMean()))));
    };

    src->traverse(process_bundle);
}

template <typename T>
inline void from(
        const typename cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>::Ptr &src,
        cslibs_ndt_3d::DistributionArray::Ptr &dst,
        const typename cslibs_gridmaps::utility::InverseModel<T>::Ptr &ivm,
        const T &threshold = 0.169)
{
    if (!src)
        return;
    src->allocatePartiallyAllocatedBundles();

    using src_map_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap<T>;
    using dst_map_t = cslibs_ndt_3d::DistributionArray;
    dst.reset(new dst_map_t());

    using point_t               = typename src_map_t::point_t;
    using distribution_t        = typename src_map_t::distribution_t;
    using distribution_bundle_t = typename src_map_t::distribution_bundle_t;
    auto sample = [&ivm](const distribution_t *d,
                         const point_t &p) -> T {
        auto evaluate = [&ivm, d, p] {
            const auto &handle = d;
            return d && handle->getDistribution() ?
                        handle->getDistribution()->sampleNonNormalized(p) * handle->getOccupancy(ivm) : 0.0;
        };
        return d ? evaluate() : 0.0;
    };
    auto sample_bundle = [&sample](const distribution_bundle_t &b,
                                   const point_t &p) -> T {
        return 0.125 * (sample(b.at(0), p) +
                        sample(b.at(1), p) +
                        sample(b.at(2), p) +
                        sample(b.at(3), p) +
                        sample(b.at(4), p) +
                        sample(b.at(5), p) +
                        sample(b.at(6), p) +
                        sample(b.at(7), p));
    };    

    using index_t = std::array<int, 3>;
    auto process_bundle = [&dst, &ivm, &threshold, &sample_bundle](const index_t &bi, const distribution_bundle_t &b) {
        typename distribution_t::distribution_t d;
        T occupancy = 0.0;

        for (std::size_t i = 0 ; i < 8 ; ++i) {
            const auto &handle = b.at(i);
            occupancy += 0.125 * handle->getOccupancy(ivm);
            if (const auto &d_tmp = handle->getDistribution())
                d += *d_tmp;
        }
        if (d.getN() == 0 || occupancy < threshold)
            return;

        dst->data.emplace_back(from(d, b.id(), sample_bundle(b, point_t(d.getMean()))));
    };
    src->traverse(process_bundle);
}
}
}

#endif // CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
