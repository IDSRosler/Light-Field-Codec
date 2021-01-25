
#include <gtest/gtest.h>
#include <numeric>
#include "Transform.h"
#include "deprecated/Transform.h"
#include "utils.h"

TEST(PartialTransformTest, zero_offset)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);

    float block[SIZE] = {0};
    float t1_block[SIZE] = {0};
    float t2_block[SIZE] = {0};

    std::iota(block, block + SIZE, 0);
    Transform tx(shape);
    old::Transform old_tx(shape);

    tx.md_forward(Transform::DCT_II, block, t1_block, {0, 0, 0, 0}, shape);
    old_tx.dct_4d(block, t2_block, shape, shape);

    ASSERT_TRUE(std::equal(t1_block, t1_block + SIZE, t2_block, [](float left, float right) { return std::round(left - right) == 0; }));
}

TEST(PartialTransformTest, smaller_shape_nicely_positioned)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    FOREACH_4D_IDX(i, x_shape, stride)
        block[i] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, {0, 0, 0, 0}, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, {0, 0, 0, 0}, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}

TEST(PartialTransformTest, smaller_shape_at_corner_x)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D x_offset(HALF_SIDE, 0, 0, 0);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    int offset = calc_offset(x_offset, stride);
    FOREACH_4D_IDX(i, x_shape, stride)
        block[i + offset] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}

TEST(PartialTransformTest, smaller_shape_at_corner_y)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D x_offset(0, HALF_SIDE, 0, 0);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    int offset = calc_offset(x_offset, stride);
    FOREACH_4D_IDX(i, x_shape, stride)
                        block[i + offset] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}

TEST(PartialTransformTest, smaller_shape_at_corner_u)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D x_offset(0, 0, HALF_SIDE, 0);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    int offset = calc_offset(x_offset, stride);
    FOREACH_4D_IDX(i, x_shape, stride)
                        block[i + offset] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}

TEST(PartialTransformTest, smaller_shape_at_corner_v)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D x_offset(0, 0, 0, HALF_SIDE);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    int offset = calc_offset(x_offset, stride);
    FOREACH_4D_IDX(i, x_shape, stride)
                        block[i + offset] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}

TEST(PartialTransformTest, trimm_bigger_shape)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t BIGGER_SIZE = SIDE + 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(BIGGER_SIZE, BIGGER_SIZE, BIGGER_SIZE, BIGGER_SIZE);
    Point4D x_offset(0, 0, 0, 0);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    std::iota(block, block + SIZE, 0);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}


TEST(PartialTransformTest, trimm_out_of_bounds_shape)
{
    constexpr std::size_t SIDE = 6;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;
    Point4D shape(SIDE, SIDE, SIDE, SIDE);
    Point4D x_shape(HALF_SIDE, HALF_SIDE, HALF_SIDE, HALF_SIDE);
    Point4D x_offset(HALF_SIDE + 1, HALF_SIDE + 1, HALF_SIDE + 1, HALF_SIDE + 1);
    Point4D correct_shape(SIDE - HALF_SIDE - 1, SIDE - HALF_SIDE - 1, SIDE - HALF_SIDE - 1, SIDE - HALF_SIDE - 1);
    Point4D stride = make_stride(shape);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    int counter = 0;
    int offset = calc_offset(x_offset, stride);
    FOREACH_4D_IDX(i, correct_shape, stride)
        block[i + offset] = static_cast<float>(counter++);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, x_offset, x_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, x_offset, x_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}


TEST(PartialTransformTest, mix_two_transforms)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;

    Point4D shape(SIDE, SIDE, SIDE, SIDE);

    Point4D dct_shape(HALF_SIDE, SIDE, SIDE, SIDE);
    Point4D dct_offset(0, 0, 0, 0);

    Point4D dst_shape(HALF_SIDE, SIDE, SIDE, SIDE);
    Point4D dst_offset(HALF_SIDE, 0, 0, 0);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    std::iota(block, block + SIZE, 0);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_block, dct_offset, dct_shape);
    tx.md_forward(Transform::DST_I, block, t_block, dst_offset, dst_shape);

    tx.md_inverse(Transform::DST_I, t_block, r_block, dst_offset, dst_shape);
    tx.md_inverse(Transform::DCT_II, t_block, r_block, dct_offset, dct_shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(block, block + SIZE, r_block));
}


TEST(PartialTransformTest, transfrom_from_descriptor)
{
    constexpr std::size_t SIDE = 8;
    constexpr std::size_t HALF_SIDE = SIDE >> 1;
    constexpr std::size_t SIZE = SIDE * SIDE * SIDE * SIDE;

    Point4D shape(SIDE, SIDE, SIDE, SIDE);

    Point4D quarter_shape(HALF_SIDE, HALF_SIDE, SIDE, SIDE);

    Point4D offset_quad_1(0, 0, 0, 0);
    Point4D offset_quad_2(HALF_SIDE, 0, 0, 0);
    Point4D offset_quad_3(0, HALF_SIDE, 0, 0);
    Point4D offset_quad_4(HALF_SIDE, HALF_SIDE, 0, 0);

    float block[SIZE] = {0};
    float t_descriptor_block[SIZE] = {0};
    float t_ground_truth_block[SIZE] = {0};

    std::iota(block, block + SIZE, 0);

    Transform tx(shape);
    tx.md_forward(Transform::DCT_II, block, t_ground_truth_block, offset_quad_1, quarter_shape);
    tx.md_forward(Transform::DST_I, block, t_ground_truth_block, offset_quad_2, quarter_shape);
    tx.md_forward(Transform::DCT_II, block, t_ground_truth_block, offset_quad_3, quarter_shape);
    tx.md_forward(Transform::DST_I, block, t_ground_truth_block, offset_quad_4, quarter_shape);

    tx.forward("s1212", block, t_descriptor_block, shape);

    ASSERT_TRUE(std::equal(t_descriptor_block, t_descriptor_block + SIZE, t_ground_truth_block));
}



TEST(PartialTransformTest, transfrom_and_reconstruct_from_descriptor)
{

    constexpr std::size_t SIZE = 15 * 15 * 13 * 13;

    Point4D shape(15, 15, 13, 13);

    float block[SIZE] = {0};
    float t_block[SIZE] = {0};
    float r_block[SIZE] = {0};

    std::iota(block, block + SIZE, 0);

    Transform tx(shape);
    const auto descriptor = "s11s11111";

    tx.forward(descriptor, block, t_block, shape);
    tx.inverse(descriptor, t_block, r_block, shape);

    std::transform(r_block, r_block + SIZE, r_block, std::roundf);
    ASSERT_TRUE(std::equal(r_block, r_block + SIZE, block));
}
