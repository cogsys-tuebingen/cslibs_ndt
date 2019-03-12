#include <gtest/gtest.h>

#include <cslibs_ndt/map/utility.hpp>
#include <cslibs_math/random/random.hpp>

const std::size_t NUM_SAMPLES = 1000;
using rng_t = cslibs_math::random::Uniform<double,1>;

inline bool get(const std::array<bool,7>& values, const std::size_t &counter)
{
    return values[counter];
}
template <std::size_t... counter>
inline bool valid(const std::array<bool,7>& values,
                  cslibs_ndt::map::detail::integer_sequence<std::size_t,counter...> counts)
{
    return cslibs_ndt::map::detail::merge<cslibs_ndt::map::detail::bool_and>(get(values,counter)...);
}

TEST(Test_cslibs_ndt, testMergeAnd)
{
    using bool_list_t = std::array<bool,7>;

    rng_t rng(-100.0, +100.0);
    for (std::size_t i=0; i<NUM_SAMPLES; ++i) {
        const bool_list_t values = {{
            rng.get() > 0,
            rng.get() > 0,
            rng.get() > 0,
            rng.get() > 0,
            rng.get() > 0,
            rng.get() > 0,
            rng.get() > 0,
        }};

        const bool ground_truth =
                values[0] &&
                values[1] &&
                values[2] &&
                values[3] &&
                values[4] &&
                values[5] &&
                values[6];

        auto is_valid = [&values]() {
            return valid(values,cslibs_ndt::map::detail::make_integer_sequence<std::size_t,std::tuple_size<bool_list_t>::value>{});
        };
        const bool merge_and_value = is_valid();
        EXPECT_EQ(merge_and_value, ground_truth);
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
