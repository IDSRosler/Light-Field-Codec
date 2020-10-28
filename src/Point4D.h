
#ifndef POINT4D_H
#define POINT4D_H

#include <iomanip>
#include <iostream>
#include <vector>

class Point4D
{
public:
  using value_type = std::size_t;
  value_type x{8}, y{8}, u{8}, v{8};
  value_type nSamples{0};

  Point4D() = default;

  Point4D(value_type x, value_type y, value_type u, value_type v);
  template <typename T>
  Point4D(const T *array);

  void updateNSamples();
  value_type getNSamples() const;

  void set(value_type x, value_type y, value_type u, value_type v);

  bool is_equal(value_type _x, value_type _y, value_type _u, value_type _v) {
      return *this == Point4D(_x, _y, _u, _v);
  }

  std::vector<value_type> to_vector() const;

  friend std::ostream &operator<<(std::ostream &os, Point4D const &point);

  friend bool operator!=(const Point4D &p1, const Point4D &p2);

  friend bool operator==(const Point4D &p1, const Point4D &p2);

  friend Point4D operator>>(const Point4D &p, const int val);

  const value_type &operator[](size_t dimension) const;

  value_type &operator[](size_t dimension);

  friend Point4D operator+(const Point4D &p, const Point4D &q);

  friend bool operator<(const Point4D &p, const Point4D &q);

  friend bool operator>(const Point4D &p, const Point4D &q);

  friend Point4D operator-(const Point4D &p, const Point4D &q);
};

struct Point4DHasher
{
  Point4D::value_type operator()(const Point4D &p) const
  {
    return p.x ^ ((p.y << 1) ^ ((p.u << 2) ^ (p.v << 3)));
  }
};

template <typename T>
Point4D::Point4D(const T *array) : x(array[0]), y(array[1]), u(array[2]), v(array[3]) { updateNSamples(); }

#endif // POINT4D_H
