#include "Point4D.h"
#include "Transform.h"
#include <algorithm>
#include <cstdlib>
#include <gtest/gtest.h>
#include <iostream>

#include <cassert>
#include <cmath>

inline int offset(int x, int y, int u, int v, Point4D &stride) {
  return x * stride.x + y * stride.y + u * stride.u + v * stride.v;
}

void flip_axis(float *block, unsigned to_flip, Point4D shape, Point4D stride) {
  auto flat_size = shape.getNSamples();
  float _block[flat_size];
  for (int v = 0; v < shape.v; v++) {
    for (int u = 0; u < shape.u; u++) {
      for (int y = 0; y < shape.y; y++) {
        for (int x = 0; x < shape.x; x++) {
          auto dx = to_flip & Transform::AXIS_X ? shape.x - 1 - 2 * x : 0;
          auto dy = to_flip & Transform::AXIS_Y ? shape.y - 1 - 2 * y : 0;
          auto du = to_flip & Transform::AXIS_U ? shape.u - 1 - 2 * u : 0;
          auto dv = to_flip & Transform::AXIS_V ? shape.v - 1 - 2 * v : 0;
          auto f_offset = offset(x, y, u, v, stride);
          auto r_offset = offset(x + dx, y + dy, u + du, v + dv, stride);
          _block[r_offset] = block[f_offset];
        }
      }
    }
  }
  for (int i = 0; i < flat_size; i++)
    block[i] = _block[i];
}

TEST(AxisFlipTest, flip_axis_x) {
  const unsigned SIZE = 16;
  float arr[SIZE];
  float rev[SIZE];
  Point4D stride(1, 0, 0, 0);
  Point4D shape(16, 1, 1, 1);
  for (int i = 0; i < SIZE; i++)
    rev[SIZE - 1 - i] = arr[i] = i;
  flip_axis(arr, Transform::AXIS_X, shape, stride);
  for (int i = 0; i < SIZE; i++)
    ASSERT_EQ(arr[i], rev[i]) << i;
}
TEST(AxisFlipTest, flip_axis_y) {
  const unsigned SIZE = 16;
  float arr[SIZE];
  float rev[SIZE];
  Point4D stride(1, 2, 4, 8);
  Point4D shape(2, 2, 2, 2);
  int counter = 0;
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto left = offset(x, y, u, v, stride);
          auto right = offset(x, shape.y - 1 - y, u, v, stride);
          rev[right] = arr[left] = counter;
          counter++;
        }
  flip_axis(arr, Transform::AXIS_Y, shape, stride);
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto index = offset(x, y, u, v, stride);
          ASSERT_EQ(rev[index], arr[index]) << x << y << u << v;
        }
}
TEST(AxisFlipTest, flip_axis_u) {
  const unsigned SIZE = 16;
  float arr[SIZE];
  float rev[SIZE];
  Point4D stride(1, 2, 4, 8);
  Point4D shape(2, 2, 2, 2);
  int counter = 0;
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto left = offset(x, y, u, v, stride);
          auto right = offset(x, y, shape.u - 1 - u, v, stride);
          rev[right] = arr[left] = counter;
          counter++;
        }
  flip_axis(arr, Transform::AXIS_U, shape, stride);
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto index = offset(x, y, u, v, stride);
          ASSERT_EQ(rev[index], arr[index]) << x << y << u << v;
        }
}
TEST(AxisFlipTest, flip_axis_v) {
  const unsigned SIZE = 16;
  float arr[SIZE];
  float rev[SIZE];
  Point4D stride(1, 2, 4, 8);
  Point4D shape(2, 2, 2, 2);
  int counter = 0;
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto left = offset(x, y, u, v, stride);
          auto right = offset(x, y, u, shape.v - 1 - v, stride);
          rev[right] = arr[left] = counter;
          counter++;
        }
  flip_axis(arr, Transform::AXIS_V, shape, stride);
  for (int v = 0; v < shape.v; v++)
    for (int u = 0; u < shape.u; u++)
      for (int y = 0; y < shape.y; y++)
        for (int x = 0; x < shape.x; x++) {
          auto index = offset(x, y, u, v, stride);
          ASSERT_EQ(rev[index], arr[index]) << x << y << u << v;
        }
}