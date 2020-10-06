#include <algorithm>
#include <cstdlib>
#include <gtest/gtest.h>

#include "Transform.h"
#include "deprecated/Transform.h"
#include "utils.h"

TEST(BackwardsCompatibilityWithOldAPI, new_dct_4d_is_consistent) {
  size_t FLAT_SIZE = 3 * 4 * 5 * 6;
  Point4D shape(3, 4, 5, 6);

  Transform new_transform(shape);
  old::Transform old_transform(shape);


  float input[FLAT_SIZE];
  float output[FLAT_SIZE];
  float expected[FLAT_SIZE];

  // Seed the random number generator
  std::srand(0);

  // Populate input with random numbers and initialize output arrays.
  for (int i = 0; i < FLAT_SIZE; i++) {
    float value = 256 * (std::rand() / (RAND_MAX * 1.0));
    input[i] = value;
    output[i] = 0;
    expected[i] = 0;
  }

  old_transform.dct_4d(input, expected, shape, shape);
  new_transform.md_forward(Transform::DCT_II, input, output, shape);

  for (int i = 0; i < FLAT_SIZE; i++)
    ASSERT_NEAR(expected[i], output[i], std::abs(expected[i] * 1e-3));
}

TEST(BackwardsCompatibilityWithOldAPI, new_idct_4d_is_consistent) {
  size_t FLAT_SIZE = 3 * 4 * 5 * 6;
  Point4D shape(3, 4, 5, 6);

  Transform new_transform(shape);
  old::Transform old_transform(shape);
  

  float input[FLAT_SIZE];
  float temp1[FLAT_SIZE];
  float temp2[FLAT_SIZE];
  float output[FLAT_SIZE];
  float expected[FLAT_SIZE];

  // Seed the random number generator
  std::srand(0);

  // Populate input with random numbers and initialize output arrays.
  for (int i = 0; i < FLAT_SIZE; i++) {
    float value = 256 * (std::rand() / (RAND_MAX * 1.0));
    input[i] = value;
    output[i] = 0;
    expected[i] = 0;
  }

  old_transform.dct_4d(input, temp1, shape, shape);
  old_transform.idct_4d(temp1, expected, shape, shape);

  new_transform.md_forward(Transform::DCT_II, input, temp2, shape);
  new_transform.md_inverse(Transform::DCT_II, temp2, output, shape);

  for (int i = 0; i < FLAT_SIZE; i++)
    ASSERT_NEAR(expected[i], output[i], std::abs(expected[i] * 1e-3));
}

TEST(BackwardsCompatibilityWithOldAPI, partial_size_dct) {
  auto MAGIC_NUMBER = 123456.0;
  size_t FLAT_SIZE = 38025;
  Point4D shape(15,15,13,13);
  Point4D smaller_shape(15,15,13,10);
  Point4D stride(1,15,225,2925);

  Transform new_transform(shape);
  old::Transform old_transform(shape);


  float input[FLAT_SIZE];
  float temp1[FLAT_SIZE];
  float temp2[FLAT_SIZE];
  float output[FLAT_SIZE];
  float expected[FLAT_SIZE];

  // Populate input with random numbers and initialize output arrays.
  for (int i = 0; i < FLAT_SIZE; i++) {
    input[i] = i;
    output[i] = expected[i] = temp1[i] = temp2[i] = MAGIC_NUMBER;
  }

  old_transform.dct_4d(input, temp1, smaller_shape, shape);
  new_transform.md_forward(Transform::DCT_II, input, temp2, smaller_shape);

  old_transform.idct_4d(temp1, expected, smaller_shape, shape);
  new_transform.md_inverse(Transform::DCT_II, temp2, output, smaller_shape);

  for (int v = 0; v < shape.v; v++) {
    for (int u = 0; u < shape.u; u++) {
      for (int y = 0; y < shape.y; y++) {
        for (int x = 0; x < shape.x; x++) {
          auto index = calc_offset(x, y, u, v, stride);
          if (!(v < smaller_shape.v && u < smaller_shape.u &&
                y < smaller_shape.y && x < smaller_shape.x)) {
            ASSERT_EQ(temp1[index], MAGIC_NUMBER);
            ASSERT_EQ(temp2[index], MAGIC_NUMBER);
            ASSERT_EQ(output[index], MAGIC_NUMBER);
            ASSERT_EQ(expected[index], MAGIC_NUMBER);
          } else {
            ASSERT_NEAR(temp1[index], temp2[index], 1);
            ASSERT_NEAR(output[index], expected[index], 1);
          }
        }
      }
    }
  }
}


TEST(BackwardsCompatibilityWithOldAPI, partial_size_dst2) {
  auto MAGIC_NUMBER = 0xabcde;
  size_t FLAT_SIZE = 38025;
  Point4D shape(15,15,13,13);
  Point4D smaller_shape(10,15,13,13);
  Point4D stride(1,15,225,2925);

  Transform new_transform(shape);


  float input[FLAT_SIZE];
  float transformed_signal[FLAT_SIZE];
  float output[FLAT_SIZE];

  // Populate input with random numbers and initialize output arrays.
  for (int i = 0; i < FLAT_SIZE; i++) {
    input[i] = i;
    output[i] = transformed_signal[i] = MAGIC_NUMBER;
  }

  new_transform.md_forward(Transform::DST_I, input, transformed_signal, smaller_shape);
  new_transform.md_inverse(Transform::DST_I, transformed_signal, output, smaller_shape);

  for (int v = 0; v < shape.v; v++) {
    for (int u = 0; u < shape.u; u++) {
      for (int y = 0; y < shape.y; y++) {
        for (int x = 0; x < shape.x; x++) {
          auto index = calc_offset(x, y, u, v, stride);
          if (!(v < smaller_shape.v && u < smaller_shape.u &&
                y < smaller_shape.y && x < smaller_shape.x)) {
            ASSERT_EQ(transformed_signal[index], MAGIC_NUMBER);
            ASSERT_EQ(output[index], MAGIC_NUMBER);
          } else {
            ASSERT_NEAR(output[index], input[index], 1);
          }
        }
      }
    }
  }
}