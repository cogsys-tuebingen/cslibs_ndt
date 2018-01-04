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

}

TEST(Test_cslibs_ndt_2d, testStaticGridmapSerialization)
{

}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
