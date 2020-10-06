
#ifndef POINT4D_H
#define POINT4D_H

//#include <zconf.h>
#include <iostream>
#include <iomanip>


class Point4D {

public:
    uint x{8}, y{8}, u{8}, v{8};
    uint nSamples{0};

    Point4D() = default;

    Point4D(uint x, uint y, uint u, uint v);

    void updateNSamples();

    uint getNSamples() const;

    friend std::ostream &operator<<(std::ostream &os, Point4D const &point);

    friend bool operator!=(const Point4D &p1, const Point4D &p2);

    friend bool operator==(const Point4D &p1, const Point4D &p2);

};

#endif //POINT4D_H
