
#include "Point4D.h"

Point4D::Point4D(uint x, uint y, uint u, uint v) : x(x), y(y), u(u), v(v) {
    this->updateNSamples();
}

Point4D::Point4D(const uint *array) 
: Point4D(array[X_AXIS], array[Y_AXIS], array[U_AXIS], array[V_AXIS])
{
}

Point4D::Point4D(const float *array) 
: Point4D(array[X_AXIS], array[Y_AXIS], array[U_AXIS], array[V_AXIS])
{
}

uint Point4D::getNSamples() const {
    return this->nSamples;
}

void Point4D::updateNSamples() {
    this->nSamples = this->x * this->y * this->u * this->v;
}

std::ostream &operator<<(std::ostream &os, Point4D const &point) {
    os << "{"
       << "x: " << std::setfill('0') << std::setw(4) << point.x << ", "
       << "y: " << std::setfill('0') << std::setw(4) << point.y << ", "
       << "u: " << std::setfill('0') << std::setw(4) << point.u << ", "
       << "v: " << std::setfill('0') << std::setw(4) << point.v
       << "}";

    return os;
}

bool operator==(const Point4D &p1, const Point4D &p2) {
    return (p1.x == p2.x && p1.y == p2.y && p1.u == p2.u && p1.v == p2.v);
}

bool operator!=(const Point4D &p1, const Point4D &p2) {
    return !(p1 == p2);
}

Point4D operator-(const Point4D &p1, const Point4D &p2) {
    if (p2.x > p1.x || p2.y > p1.y || p2.u > p1.u || p2.v > p1.v) 
        throw std::overflow_error("Tried to subtract a bigger value from "
                                   "smaller one.");

    Point4D p(p1.x - p2.x, p1.y - p2.y, p1.u - p2.u, p1.v - p2.v);
    return p;
}

Point4D operator+(const Point4D &p1, const Point4D &p2) {
    Point4D p(p1.x + p2.x, p1.y + p2.y, p1.u + p2.u, p1.v + p2.v);
    return p;
}


Point4D operator+(const Point4D &p, const uint delta) {
    Point4D q(delta, delta, delta, delta);
    return p + q;
}

Point4D operator-(const Point4D &p, const uint delta) {
    Point4D q(delta, delta, delta, delta);
    return p - q;
}

size_t *Point4D::to_array() const {
    size_t *array = new size_t[4];
    array[X_AXIS] = x;
    array[Y_AXIS] = y;
    array[U_AXIS] = u;
    array[V_AXIS] = v;

    return array;
}