
#include "Point4D.h"

Point4D::Point4D(uint x, uint y, uint u, uint v) : x(x), y(y), u(u), v(v) {
    this->updateNSamples();
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



