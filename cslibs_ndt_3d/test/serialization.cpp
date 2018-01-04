#include <gtest/gtest.h>

#include <cslibs_ndt_3d/serialization/dynamic_maps/gridmap.hpp>
#include <cslibs_ndt_3d/serialization/static_maps/gridmap.hpp>

#include <cslibs_math/random/random.hpp>

const std::size_t MIN_NUM_SAMPLES = 100;
const std::size_t MAX_NUM_SAMPLES = 1000;

template <std::size_t Dim>
using rng_t = typename cslibs_math::random::Uniform<Dim>;

TEST(Test_cslibs_ndt_3d, testDynamicGridmapSerialization)
{
    using map_t = cslibs_ndt_3d::dynamic_maps::Gridmap;
    rng_t<1> rng_coord(-10.0, 10.0);
    rng_t<1> rng_angle(-M_PI, M_PI);

    // fill map
    cslibs_math_3d::Transform3d origin(cslibs_math_3d::Vector3d(rng_coord.get(), rng_coord.get(), rng_coord.get()),
                                       cslibs_math_3d::Quaternion(rng_angle.get(), rng_angle.get(), rng_angle.get()));
    const double resolution = rng_t<1>(1.0, 5.0).get();
    typename map_t::Ptr map(new map_t(origin, resolution));
    const int num_samples = static_cast<int>(rng_t<1>(MIN_NUM_SAMPLES, MAX_NUM_SAMPLES).get());
    for (int i = 0 ; i < num_samples ; ++ i) {
        const cslibs_math_3d::Point3d p(rng_coord.get(), rng_coord.get(), rng_coord.get());
        map->add(p);
    }

    // serialization
    YAML::Node n(map);

    // de-serialization
    const typename map_t::Ptr & map_converted = n.as<typename map_t::Ptr>();
    EXPECT_NE(map_converted, nullptr);

    // tests
    EXPECT_NEAR(map->getResolution(),       map_converted->getResolution(),       1e-9);
    EXPECT_NEAR(map->getBundleResolution(), map_converted->getBundleResolution(), 1e-9);
    EXPECT_NEAR(map->getHeight(),           map_converted->getHeight(),           1e-9);
    EXPECT_NEAR(map->getWidth(),            map_converted->getWidth(),            1e-9);

    EXPECT_EQ(map->getMinDistributionIndex()[0], map_converted->getMinDistributionIndex()[0]);
    EXPECT_EQ(map->getMinDistributionIndex()[1], map_converted->getMinDistributionIndex()[1]);
    EXPECT_EQ(map->getMinDistributionIndex()[2], map_converted->getMinDistributionIndex()[2]);
    EXPECT_EQ(map->getMaxDistributionIndex()[0], map_converted->getMaxDistributionIndex()[0]);
    EXPECT_EQ(map->getMaxDistributionIndex()[1], map_converted->getMaxDistributionIndex()[1]);
    EXPECT_EQ(map->getMaxDistributionIndex()[2], map_converted->getMaxDistributionIndex()[2]);

    EXPECT_NEAR(map->getMin()(0), map_converted->getMin()(0), 1e-9);
    EXPECT_NEAR(map->getMin()(1), map_converted->getMin()(1), 1e-9);
    EXPECT_NEAR(map->getMin()(2), map_converted->getMin()(2), 1e-9);
    EXPECT_NEAR(map->getMax()(0), map_converted->getMax()(0), 1e-9);
    EXPECT_NEAR(map->getMax()(1), map_converted->getMax()(1), 1e-9);
    EXPECT_NEAR(map->getMax()(2), map_converted->getMax()(2), 1e-9);

    EXPECT_NEAR(map->getOrigin().tx(),           map_converted->getOrigin().tx(),           1e-9);
    EXPECT_NEAR(map->getOrigin().ty(),           map_converted->getOrigin().ty(),           1e-9);
    EXPECT_NEAR(map->getOrigin().tz(),           map_converted->getOrigin().tz(),           1e-9);
    EXPECT_NEAR(map->getOrigin().roll(),         map_converted->getOrigin().roll(),         1e-9);
    EXPECT_NEAR(map->getOrigin().pitch(),        map_converted->getOrigin().pitch(),        1e-9);
    EXPECT_NEAR(map->getOrigin().yaw(),          map_converted->getOrigin().yaw(),          1e-9);
    EXPECT_NEAR(map->getInitialOrigin().tx(),    map_converted->getInitialOrigin().tx(),    1e-9);
    EXPECT_NEAR(map->getInitialOrigin().ty(),    map_converted->getInitialOrigin().ty(),    1e-9);
    EXPECT_NEAR(map->getInitialOrigin().tz(),    map_converted->getInitialOrigin().tz(),    1e-9);
    EXPECT_NEAR(map->getInitialOrigin().roll(),  map_converted->getInitialOrigin().roll(),  1e-9);
    EXPECT_NEAR(map->getInitialOrigin().pitch(), map_converted->getInitialOrigin().pitch(), 1e-9);
    EXPECT_NEAR(map->getInitialOrigin().yaw(),   map_converted->getInitialOrigin().yaw(),   1e-9);

    for (int idx = map->getMinDistributionIndex()[0] ; idx <= map->getMaxDistributionIndex()[0] ; ++ idx) {
        for (int idy = map->getMinDistributionIndex()[1] ; idy <= map->getMaxDistributionIndex()[1] ; ++ idy) {
            for (int idz = map->getMinDistributionIndex()[2] ; idz <= map->getMinDistributionIndex()[2] ; ++ idz) {
                std::array<int, 3> bi({idx, idy, idz});
                if (typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_bundle_t* b =
                        map->getDistributionBundle(bi)) {
                    typename cslibs_ndt_3d::dynamic_maps::Gridmap::distribution_bundle_t* bb =
                            map_converted->getDistributionBundle(bi);
                    EXPECT_NE(bb, nullptr);

                    for (std::size_t i = 0 ; i < 8 ; ++ i) {
                        EXPECT_NE(b->at(i),  nullptr);
                        EXPECT_NE(bb->at(i), nullptr);

                        const cslibs_math::statistics::Distribution<3, 3> & d  = b->at(i)->getHandle()->data();
                        const cslibs_math::statistics::Distribution<3, 3> & dd = bb->at(i)->getHandle()->data();
                        EXPECT_EQ(d.getN(), dd.getN());

                        for (std::size_t j = 0 ; j < 3 ; ++ j) {
                            EXPECT_NEAR(d.getMean()(j), dd.getMean()(j), 1e-3);
                            for (std::size_t k = 0 ; k < 3 ; ++ k) {
                                EXPECT_NEAR(d.getCorrelated()(j, k), dd.getCorrelated()(j, k), 1e-3);
                                EXPECT_NEAR(d.getCovariance()(j, k), dd.getCovariance()(j, k), 1e-3);
                                EXPECT_NEAR(d.getInformationMatrix()(j, k), dd.getInformationMatrix()(j, k), 1e-3);
                            }
                        }
                    }
                } else
                    EXPECT_EQ(map_converted->getDistributionBundle(bi), nullptr);
            }
        }
    }
}

TEST(Test_cslibs_ndt_3d, testStaticGridmapSerialization)
{/*
    using map_t = cslibs_ndt_3d::static_maps::Gridmap;
    rng_t<1> rng_coord(-10.0, 10.0);
    rng_t<1> rng_angle(-M_PI, M_PI);
    rng_t<1> rng_size(100.0, 200.0);
/*
    // fill map
    cslibs_math_3d::Transform3d origin(rng_coord.get(), rng_coord.get(), rng_coord.get(),
                                       rng_angle.get(), rng_angle.get(), rng_angle.get());
    const double resolution = rng_t<1>(1.0, 5.0).get();
    const std::size_t size_x = rng_size.get();
    const std::size_t size_y = rng_size.get();
    typename map_t::Ptr map(new map_t(origin, resolution, {{size_x, size_y}}));
    const int num_samples = static_cast<int>(rng_t<1>(MIN_NUM_SAMPLES, MAX_NUM_SAMPLES).get());

    const double max_coord_x = static_cast<double>(size_x) * resolution;
    const double max_coord_y = static_cast<double>(size_y) * resolution;
    rng_t<1> rng_coord_x(0.0, max_coord_x);
    rng_t<1> rng_coord_y(0.0, max_coord_y);
    for (int i = 0 ; i < num_samples ; ++ i) {
        const cslibs_math_2d::Point2d p(rng_coord_x.get(), rng_coord_y.get());
        map->add(origin * p);
    }

    // serialization
    YAML::Node n(map);

    // de-serialization
    const typename map_t::Ptr & map_converted = n.as<typename map_t::Ptr>();
    EXPECT_NE(map_converted, nullptr);

    // tests
    EXPECT_NEAR(map->getResolution(),       map_converted->getResolution(),       1e-9);
    EXPECT_NEAR(map->getBundleResolution(), map_converted->getBundleResolution(), 1e-9);
    EXPECT_NEAR(map->getHeight(),           map_converted->getHeight(),           1e-9);
    EXPECT_NEAR(map->getWidth(),            map_converted->getWidth(),            1e-9);

    EXPECT_EQ(map->getSize()[0],       map_converted->getSize()[0]);
    EXPECT_EQ(map->getSize()[1],       map_converted->getSize()[1]);
    EXPECT_EQ(map->getBundleSize()[0], map_converted->getBundleSize()[0]);
    EXPECT_EQ(map->getBundleSize()[1], map_converted->getBundleSize()[1]);

    EXPECT_NEAR(map->getOrigin().tx(),  map_converted->getOrigin().tx(),  1e-9);
    EXPECT_NEAR(map->getOrigin().ty(),  map_converted->getOrigin().ty(),  1e-9);
    EXPECT_NEAR(map->getOrigin().yaw(), map_converted->getOrigin().yaw(), 1e-9);

    for (int idx = 0 ; idx < static_cast<int>(map->getBundleSize()[0]) ; ++ idx) {
        for (int idy = 0 ; idy < static_cast<int>(map->getBundleSize()[1]) ; ++ idy) {
            std::array<int, 2> bi({idx, idy});
            if (typename cslibs_ndt_2d::static_maps::Gridmap::distribution_bundle_t* b =
                    map->getDistributionBundle(bi)) {
                typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_bundle_t* bb =
                        map_converted->getDistributionBundle(bi);
                EXPECT_NE(bb, nullptr);

                for (std::size_t i = 0 ; i < 4 ; ++ i) {
                    EXPECT_NE(b->at(i),  nullptr);
                    EXPECT_NE(bb->at(i), nullptr);

                    const cslibs_math::statistics::Distribution<2, 3> & d  = b->at(i)->getHandle()->data();
                    const cslibs_math::statistics::Distribution<2, 3> & dd = bb->at(i)->getHandle()->data();
                    EXPECT_EQ(d.getN(), dd.getN());

                    for (std::size_t j = 0 ; j < 2 ; ++ j) {
                        EXPECT_NEAR(d.getMean()(j), dd.getMean()(j), 1e-3);
                        for (std::size_t k = 0 ; k < 2 ; ++ k) {
                            EXPECT_NEAR(d.getCorrelated()(j, k), dd.getCorrelated()(j, k), 1e-3);
                            EXPECT_NEAR(d.getCovariance()(j, k), dd.getCovariance()(j, k), 1e-3);
                            EXPECT_NEAR(d.getInformationMatrix()(j, k), dd.getInformationMatrix()(j, k), 1e-3);
                        }
                    }
                }
            } else
                EXPECT_EQ(map_converted->getDistributionBundle(bi), nullptr);
        }
    }*/
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
