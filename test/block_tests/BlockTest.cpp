#include "Block.h"
#include <gtest/gtest.h>
#include <stdexcept>

TEST(BlockTests, block_uses_same_memory_space_as_original_array) {
  int array[100];
  Block<int> b(array, 100);

  for (int i = 0; i < 100; i++) {
    array[i] = i;
  }

  bool same_values = true;
  for (int i = 0; i < 100; i++) {
    if (b(i) != array[i])
      same_values = false;
  }
  ASSERT_TRUE(same_values);
}

TEST(BlockTests, throw_exception_if_tries_to_access_invalid_location) {
  int variable;
  Block<int> b(&variable, 1);

  ASSERT_THROW((b(2) = 4), std::out_of_range);
}

TEST(BlockTests, reshape_block) {
  int array[100];
  Block<int> b(array, 100);

  auto shape = b.shape();
  ASSERT_EQ(shape[0], 100);
  ASSERT_EQ(shape.size(), 1);

  b.reshape({10, 10});
  auto new_shape = b.shape();

  ASSERT_EQ(new_shape[0], 10);
  ASSERT_EQ(new_shape[1], 10);
  ASSERT_EQ(new_shape.size(), 2);
}

TEST(BlockTests, index_multidimensional_block) {
  int array[100];
  Block<int> b(array, 100);

  // Initialize array
  for (int i = 0; i < 100; i++)
    array[i] = i;

  b.reshape({10, 10});
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
      ASSERT_EQ(b(i, j), array[i * 10 + j]);
}

TEST(BlockTests, reshape_invalid_size_throws_exception) {
  int var;
  Block<int> b(&var, 15);

  EXPECT_NO_THROW(b.reshape({3, 5}));
  EXPECT_NO_THROW(b.reshape({5, 3}));
  ASSERT_THROW(b.reshape({1, 10}), std::length_error);
}

TEST(BlockTests, access_block_elements_through_view) {
  using Range = Block<int>::Range;
  int array[1000];
  Block<int> b(array, 1000);
  b.reshape({10, 10, 10});

  // Initialize array
  for (int i = 0; i < 1000; i++) {
    array[i] = i;
  }

  auto v1 = b.view({Range(0, 5), Range(0, 5), Range(0, 5)}); // First Octant
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 5; k++)
        ASSERT_EQ(v1(i, j, k), array[i * 100 + j * 10 + k]);

  auto v2 = b.view({Range(6, 10), Range(0, 5), Range(0, 5)}); // Second Octant
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 5; k++)
        ASSERT_EQ(v2(i, j, k), array[(i + 6) * 100 + j * 10 + k]);

  auto v3 = b.view({Range(2, 8), Range(1, 10), Range(7, 9)}); // view 6 x 9 x 2
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 9; j++)
      for (int k = 0; k < 2; k++)
        ASSERT_EQ(v3(i, j, k), array[(i + 2) * 100 + (j + 1) * 10 + (k + 7)]);
}

TEST(BlockTests, view_with_invalid_ranges_throws_exception) {
  using Range = Block<int>::Range;
  int array[1000];
  Block<int> b(array, 1000);
  b.reshape({10, 10, 10});

  ASSERT_THROW(b.view({Range(5, 11), Range(2, 3), Range(3, 4)}),
               std::out_of_range);
  ASSERT_THROW(b.view({Range(2, 3), Range(5, 11), Range(3, 4)}),
               std::out_of_range);
  ASSERT_THROW(b.view({Range(2, 3), Range(3, 4), Range(5, 11)}),
               std::out_of_range);
  ASSERT_THROW(b.view({Range(2, 3), Range(3, 4)}), std::invalid_argument);
  ASSERT_THROW(b.view({Range(4, 0), Range(2, 3), Range(3, 4)}),
               std::invalid_argument);
}

TEST(BlockTests, create_view_of_view) {
  using Range = Block<int>::Range;
  int array[8 * 8 * 8];
  size_t x0, x1, dx, y0, y1, dy, z0, z1, dz;
  Block<int> b(array, 8 * 8 * 8);

  b.reshape({8, 8, 8});
  x0 = 3;
  x1 = 6;
  y0 = 2;
  y1 = 7;
  z0 = 1;
  z1 = 5;
  dy = y1 - y0;
  dx = x1 - x0;
  dz = z1 - z0;

  // Initialize array
  for (int i = 0; i < 8 * 8 * 8; i++)
    array[i] = i;

  // Slice b three times
  auto v1 = b.view({Range(x0, x1), Range(0, 8), Range(0, 8)});
  auto v2 = v1.view({Range(0, dx), Range(y0, y1), Range(0, 8)});
  auto v3 = v2.view({Range(0, dx), Range(0, dy), Range(z0, z1)});

  // Slice b at once.
  auto v = b.view({Range(x0, x1), Range(y0, y1), Range(z0, z1)});

  for (size_t x = 0; x < dx; x++)
    for (size_t y = 0; y < dy; y++)
      for (size_t z = 0; z < dz; z++)
        ASSERT_EQ(v3(x, y, z), v(x, y, z));
}

TEST(BlockTest, flat_size_is_constant_regardless_of_reshaping) {
  int array[120];
  Block<int> b(array, 120);

  b.reshape({3, 40});
  ASSERT_EQ(b.flat_size(), 120);
  b.reshape({40, 3});
  ASSERT_EQ(b.flat_size(), 120);
  b.reshape({5, 4, 6});
  ASSERT_EQ(b.flat_size(), 120);
  b.reshape({120});
  ASSERT_EQ(b.flat_size(), 120);
  b.reshape({1, 1, 120});
  ASSERT_EQ(b.flat_size(), 120);
}

TEST(BlockTest, forward_iterator_is_consistent) {
  int array[120];
  Block<int> b(array, 120);

  b.reshape({2, 5, 3, 4});
  int i = 0;
  for (auto it = b.begin(); it != b.end(); i++, it++) {
    *it = i * i;
    ASSERT_EQ(array[i], i * i);
  }
  ASSERT_EQ(i, 120);
}

TEST(BlockTest, reverse_iterator_is_consistent) {
  int array[120];
  Block<int> b(array, 120);

  b.reshape({2, 5, 3, 4});
  int i = 119;
  for (auto it = b.rbegin(); it != b.rend(); i--, it++) {
    *it = i * i;
    ASSERT_EQ(array[i], i * i);
  }
  ASSERT_EQ(i, -1);
}

TEST(BlockTests, reshape_block_parameter_pack) {
  int array[100];
  Block<int> b(array, 100);

  auto shape = b.shape();
  ASSERT_EQ(shape[0], 100);
  ASSERT_EQ(shape.size(), 1);

  b.reshape(10, 10);
  auto new_shape = b.shape();

  ASSERT_EQ(new_shape[0], 10);
  ASSERT_EQ(new_shape[1], 10);
  ASSERT_EQ(new_shape.size(), 2);
}

TEST(BlockTests, invalid_reshape_block_parameter_pack) {
  int array[100];
  Block<int> b(array, 100);

  ASSERT_THROW(b.reshape(5, 5), std::length_error);
}

TEST(BlockTests, cant_reshape_view) {
  using Range = Block<int>::Range;
  int array[8 * 8 * 8];
  Block b(array, 8 * 8 * 8);

  b.reshape({8, 8, 8});

  auto v = b.view({Range(0, 4), Range(0, 4), Range(0, 4)});
  ASSERT_THROW(v.reshape(8, 8), std::logic_error);
}

TEST(BlockTests, reshape_with_fake_dimension) {
  using Range = Block<int>::Range;
  int array[10];
  Block b(array, 10);

  // Initialize array
  for (int i = 0; i < 10; i++)
    array[i] = i + 1;

  b.reshape(1, 10, 1);
  ASSERT_EQ(array[5], b(0, 5, 0));
}

TEST(BlockTests, shape_of_view_is_consistent) {
  using Range = Block<int>::Range;
  int array[8 * 8 * 8];
  Block b(array, 8 * 8 * 8);

  b.reshape({8, 8, 8});

  auto v = b.view({Range(0, 4), Range(1, 4), Range(2, 7)});
  auto shape = v.shape();

  ASSERT_EQ(shape.size(), 3);
  ASSERT_EQ(shape[0], 4);
  ASSERT_EQ(shape[1], 3);
  ASSERT_EQ(shape[2], 5);
}


