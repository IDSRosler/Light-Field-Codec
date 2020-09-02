
#include "Point4D.h"
#include "sstream"

Point4D::Point4D(Point4D::value_type x, Point4D::value_type y, Point4D::value_type u, Point4D::value_type v) : x(x), y(y), u(u), v(v)
{
  this->updateNSamples();
}

Point4D::value_type Point4D::getNSamples() const { return this->nSamples; }

void Point4D::updateNSamples()
{
  this->nSamples = this->x * this->y * this->u * this->v;
}

std::ostream &operator<<(std::ostream &os, Point4D const &point)
{
  std::stringstream ss;
  ss << "("
     << point.x << ", "
     << point.y << ", "
     << point.u << ", "
     << point.v << ")";
  os << ss.str();
  return os;
}

bool operator==(const Point4D &p1, const Point4D &p2)
{
  return (p1.x == p2.x && p1.y == p2.y && p1.u == p2.u && p1.v == p2.v);
}

bool operator!=(const Point4D &p1, const Point4D &p2) { return !(p1 == p2); }

Point4D operator-(const Point4D &p1, const Point4D &p2)
{
  if (p2.x > p1.x || p2.y > p1.y || p2.u > p1.u || p2.v > p1.v)
    throw std::overflow_error("Tried to subtract a bigger value from "
                              "smaller one.");

  Point4D p(p1.x - p2.x, p1.y - p2.y, p1.u - p2.u, p1.v - p2.v);
  return p;
}

Point4D operator+(const Point4D &p1, const Point4D &p2)
{
  Point4D p(p1.x + p2.x, p1.y + p2.y, p1.u + p2.u, p1.v + p2.v);
  return p;
}

Point4D operator+(const Point4D &p, const Point4D::value_type delta)
{
  Point4D q(delta, delta, delta, delta);
  return p + q;
}

Point4D operator-(const Point4D &p, const Point4D::value_type delta)
{
  Point4D q(delta, delta, delta, delta);
  return p - q;
}

void Point4D::set(Point4D::value_type x, Point4D::value_type y, Point4D::value_type u, Point4D::value_type v)
{
  this->x = x;
  this->y = y;
  this->u = u;
  this->v = v;
  updateNSamples();
}

std::vector<Point4D::value_type> Point4D::to_vector() const
{
  std::vector<Point4D::value_type> vec;
  vec.resize(4);
  vec[0] = x;
  vec[1] = y;
  vec[2] = u;
  vec[3] = v;
  return vec;
}

Point4D operator>>(const Point4D &p, const int val)
{
  Point4D q;
  q.x = p.x >> val;
  q.y = p.y >> val;
  q.u = p.u >> val;
  q.v = p.v >> val;
  return q;
}

const Point4D::value_type &Point4D::operator[](size_t dimension) const
{
  switch (dimension)
  {
  case 0:
    return x;
  case 1:
    return y;
  case 2:
    return u;
  case 3:
    return v;
  }
  throw;
}

Point4D::value_type &Point4D::operator[](size_t dimension)
{
  switch (dimension)
  {
  case 0:
    return x;
  case 1:
    return y;
  case 2:
    return u;
  case 3:
    return v;
  }
  throw;
}

bool operator<(const Point4D &p, const Point4D &q)
{
  return p.x < q.x && p.y < q.y && p.u < q.u && p.v < q.v;
}

bool operator>(const Point4D &p, const Point4D &q)
{
  return p.x > q.x && p.y > q.y && p.u > q.u && p.v > q.v;
}