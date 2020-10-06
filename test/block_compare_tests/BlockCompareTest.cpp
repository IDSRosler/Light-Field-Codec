#include <gtest/gtest.h>
#include <random>
#include "BlockCompare.h"

struct BlockCompareTest : public ::testing::Test {   

};


TEST_F(BlockCompareTest, block_with_lower_average_wins)
{
    std::default_random_engine generator;
    std::normal_distribution<float> d1(100,5);
    std::normal_distribution<float> d2(5,5);
    float a[1000];
    float b[1000];

    for (int i = 0; i < 1000; i++) {
        a[i] = (float) d1(generator);
        b[i] = (float) d2(generator);
    }

    ASSERT_TRUE(BlockComparator<float>::compare(a, b, 1000) > 0);
}

TEST_F(BlockCompareTest, block_with_lower_coefficients_wins)
{
    float a[1000];
    float b[1000];

    for (int i = 0; i < 1000; i++) {
        a[i] = 2;
        b[i] = 10;
    }

    ASSERT_TRUE(BlockComparator<float>::compare(a, b, 1000) < 0);
}

TEST_F(BlockCompareTest, block_with_least_energy_wins)
{
    float a[1000];
    float b[1000];

    for (int i = 0; i < 1000; i++) {
        a[i] = b[i] = 1;
        if (i % 5 == 0)
            a[i] = 2;   // More multiples of 5 means more 2s
        if (i % 7 == 0)
            b[i] = 2;   
    }

    ASSERT_TRUE(BlockComparator<float>::compare(a, b, 1000) > 0);
}
