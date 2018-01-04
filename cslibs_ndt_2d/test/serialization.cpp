#include <gtest/gtest.h>

#include <cslibs_ndt_2d/serialization/dynamic_maps/gridmap.hpp>
//#include <cslibs_ndt_2d/serialization/static_maps/gridmap.hpp>

#include <cslibs_math/random/random.hpp>

const std::size_t MIN_NUM_SAMPLES = 100;
const std::size_t MAX_NUM_SAMPLES = 1000;

template <std::size_t Dim>
using rng_t = typename cslibs_math::random::Uniform<Dim>;

TEST(Test_cslibs_ndt_2d, testDynamicGridmapSerialization)
{
    using map_t = cslibs_ndt_2d::dynamic_maps::Gridmap;
    rng_t<1> rng_coord(-10.0, 10.0);

    // fill map
    cslibs_math_2d::Transform2d origin(rng_coord.get(), rng_coord.get(), rng_t<1>(-M_PI, M_PI).get());
    const double resolution = rng_t<1>(1.0, 5.0).get();
    typename map_t::Ptr map(new map_t(origin, resolution));
    const int num_samples = static_cast<int>(rng_t<1>(MIN_NUM_SAMPLES, MAX_NUM_SAMPLES).get());
    for (int i = 0 ; i < num_samples ; ++ i) {
        const cslibs_math_2d::Point2d p(rng_coord.get(), rng_coord.get());
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
    EXPECT_EQ(map->getMaxDistributionIndex()[0], map_converted->getMaxDistributionIndex()[0]);
    EXPECT_EQ(map->getMaxDistributionIndex()[1], map_converted->getMaxDistributionIndex()[1]);

    for (int idx = map->getMinDistributionIndex()[0] ; idx <= map->getMaxDistributionIndex()[0] ; ++ idx) {
        for (int idy = map->getMinDistributionIndex()[1] ; idy <= map->getMaxDistributionIndex()[1] ; ++ idy) {
            std::array<int, 2> bi({idx, idy});
            if (typename cslibs_ndt_2d::dynamic_maps::Gridmap::distribution_bundle_t* b =
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
    }
}

TEST(Test_cslibs_ndt_2d, testStaticGridmapSerialization)
{

}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
