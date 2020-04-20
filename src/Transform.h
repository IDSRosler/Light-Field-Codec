
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <zconf.h>
#include "Typedef.h"
#include "Point4D.h"
#include "TransformContext.h"
#include <cmath>

class Transform
{
private:
    Point4D stride_lightfield;
    Point4D size_lightfield;
    TransformContext<float> *ctx;

public:
    explicit Transform(const Point4D &dimBlock);
    ~Transform();

    void dct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize);
    void idct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize);
};

#endif //TRANSFORM_H
