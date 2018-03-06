#ifndef CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
#define CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP

#include <cslibs_ndt_3d/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/dynamic_maps/occupancy_gridmap.hpp>

#include <cslibs_ndt_3d/DistributionArray.h>

namespace cslibs_ndt_3d {
namespace conversion {
inline void from(
        const cslibs_ndt_3d::dynamic_maps::Gridmap::Ptr &src,
        cslibs_ndt_3d::DistributionArray::Ptr           &dst)
{
    if (!src)
        return;

    using index_t   = std::array<int, 3>;
    using point_t   = cslibs_math_3d::Point3d;
    using dst_map_t = cslibs_ndt_3d::DistributionArray;
    dst.reset(new dst_map_t());

    std::vector<index_t> indices;
    src->getBundleIndices(indices);

    using distribution_t = cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_t;
    using distribution_bundle_t = cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_bundle_t;
    auto sample = [](const distribution_t *d,
                     const point_t &p) {
        return d ? d->data().sampleNonNormalized(p) : 0.0;
    };
    auto sample_bundle = [&sample](const distribution_bundle_t *b,
                                   const point_t &p) {
        return 0.125 * (sample(b->at(0), p) +
                        sample(b->at(1), p) +
                        sample(b->at(2), p) +
                        sample(b->at(3), p) +
                        sample(b->at(4), p) +
                        sample(b->at(5), p) +
                        sample(b->at(6), p) +
                        sample(b->at(7), p));
    };

    for(auto &bi : indices)
        if (const distribution_bundle_t *b = src->getDistributionBundle(bi)) {
            distribution_t::distribution_t d;
            for (std::size_t i = 0; i < 8; ++ i)
                d += b->at(i)->getHandle()->data();

            if (d.getN() == 0)
                continue;
            point_t p(d.getMean());

            Distribution distr;
            distr.id.data = b->id();
            for (int i = 0; i < 3; ++ i) {
                distr.mean[i].data          = p(i);
                distr.eigen_values[i].data  = d.getEigenValues()(i);
            }
            for (int i = 0; i < 9; ++ i) {
                distr.eigen_vectors[i].data = d.getEigenVectors()(i);
                distr.covariance[i].data    = d.getCovariance()(i);
            }
            distr.prob.data = static_cast<double>(sample_bundle(b, p));
            dst->data.emplace_back(distr);
        }

    dst->header.stamp = ros::Time::now();
}

inline void from(
        const cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::Ptr &src,
        cslibs_ndt_3d::DistributionArray::Ptr                    &dst,
        const cslibs_gridmaps::utility::InverseModel::Ptr        &ivm)
{
    if (!src)
        return;

    using index_t   = std::array<int, 3>;
    using point_t   = cslibs_math_3d::Point3d;
    using dst_map_t = cslibs_ndt_3d::DistributionArray;
    dst.reset(new dst_map_t());

    std::vector<index_t> indices;
    src->getBundleIndices(indices);

    using distribution_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_t;
    using distribution_bundle_t = cslibs_ndt_3d::dynamic_maps::OccupancyGridmap::distribution_bundle_t;
    auto sample = [&ivm](const distribution_t *d,
                         const point_t &p) {
        return (d && d->getDistribution()) ?
                    (d->getDistribution()->sampleNonNormalized(p) * d->getOccupancy(ivm)) : 0.0;

    };
    auto sample_bundle = [&sample](const distribution_bundle_t *b,
                                   const point_t &p) {
        return 0.125 * (sample(b->at(0), p) +
                        sample(b->at(1), p) +
                        sample(b->at(2), p) +
                        sample(b->at(3), p) +
                        sample(b->at(4), p) +
                        sample(b->at(5), p) +
                        sample(b->at(6), p) +
                        sample(b->at(7), p));
    };

    for(auto &bi : indices)
        if (const distribution_bundle_t *b = src->getDistributionBundle(bi)) {
            distribution_t::distribution_t d;
            for (std::size_t i = 0; i < 8; ++ i)
                if (b->at(i)->getDistribution())
                    d += *(b->at(i)->getDistribution());

            if (d.getN() == 0)
                continue;
            point_t p(d.getMean());

            Distribution distr;
            distr.id.data = b->id();
            for (int i = 0; i < 3; ++ i) {
                distr.mean[i].data          = p(i);
                distr.eigen_values[i].data  = d.getEigenValues()(i);
            }
            for (int i = 0; i < 9; ++ i) {
                distr.eigen_vectors[i].data = d.getEigenVectors()(i);
                distr.covariance[i].data    = d.getCovariance()(i);
            }
            distr.prob.data = static_cast<double>(sample_bundle(b, p));
            dst->data.emplace_back(distr);
        }

    dst->header.stamp = ros::Time::now();
}
}
}

#endif // CSLIBS_NDT_3D_CONVERSION_DISTRIBUTIONS_HPP
