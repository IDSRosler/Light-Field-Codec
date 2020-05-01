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

  // Initialize array
  for (int i = 0; i < 100; i++) {
    array[i] = i;
  }

  b.reshape({10, 10});

  bool same_values = true;
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
      if (b(i, j) != array[i * 10 + j])
        same_values = false;

  ASSERT_TRUE(same_values);
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
  x0 = 3;  x1 = 6;
  y0 = 2;  y1 = 7;
  z0 = 1;  z1 = 5;
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