#include <gtest/gtest.h>

#include <cslibs_ndt/map/utility.hpp>
#include <cslibs_math/random/random.hpp>

const std::size_t NUM_SAMPLES = 1000;
using rng_t = cslibs_math::random::Uniform<double,1>;

TEST(Test_cslibs_ndt, testGenerateIndices2d)
{
    using index_t = std::array<int,2>;
    using index_list_t = std::array<index_t,4>;

    rng_t rng(-100.0, +100.0);
    for (std::size_t i=0; i<NUM_SAMPLES; ++i) {
        const index_t bi = {static_cast<int>(rng.get()), static_cast<int>(rng.get())};

        const int divx = cslibs_math::common::div(bi[0], 2);
        const int divy = cslibs_math::common::div(bi[1], 2);
        const int modx = cslibs_math::common::mod(bi[0], 2);
        const int mody = cslibs_math::common::mod(bi[1], 2);

        const index_list_t ground_truth_list = {{
            {divx,         divy},
            {divx + modx,  divy},
            {divx,         divy + mody},
            {divx + modx,  divy + mody}
        }};

        const index_list_t generated_list = cslibs_ndt::map::detail::generate_indices<index_list_t,2>(bi);
        for (std::size_t i=0; i<4; ++i)
            EXPECT_EQ(ground_truth_list[i], generated_list[i]);
    }
}

TEST(Test_cslibs_ndt, testGenerateIndices3d)
{
    using index_t = std::array<int,3>;
    using index_list_t = std::array<index_t,8>;

    rng_t rng(-100.0, +100.0);
    for (std::size_t i=0; i<NUM_SAMPLES; ++i) {
        const index_t bi = {static_cast<int>(rng.get()), static_cast<int>(rng.get())};

        const int divx = cslibs_math::common::div(bi[0], 2);
        const int divy = cslibs_math::common::div(bi[1], 2);
        const int divz = cslibs_math::common::div(bi[2], 2);
        const int modx = cslibs_math::common::mod(bi[0], 2);
        const int mody = cslibs_math::common::mod(bi[1], 2);
        const int modz = cslibs_math::common::mod(bi[2], 2);

        const index_list_t ground_truth_list = {{
            {divx,         divy,           divz},
            {divx + modx,  divy,           divz},
            {divx,         divy + mody,    divz},
            {divx + modx,  divy + mody,    divz},
            {divx,         divy,           divz + modz},
            {divx + modx,  divy,           divz + modz},
            {divx,         divy + mody,    divz + modz},
            {divx + modx,  divy + mody,    divz + modz}
        }};

        const index_list_t generated_list = cslibs_ndt::map::detail::generate_indices<index_list_t,3>(bi);
        for (std::size_t i=0; i<8; ++i)
            EXPECT_EQ(ground_truth_list[i], generated_list[i]);
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
