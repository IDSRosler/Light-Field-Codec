
#include <cassert>
#include "Transform.h"
#include "DiscreteCosineTransformContext4D.h"

Transform::Transform(const Point4D &dimBlock)
{
    size_lightfield = dimBlock;
    stride_lightfield.x = 1;
    stride_lightfield.y = dimBlock.x;
    stride_lightfield.u = stride_lightfield.y * dimBlock.y;
    stride_lightfield.v = stride_lightfield.u * dimBlock.u;
    ctx = new DiscreteCosineTransformContext4D<float>(size_lightfield,
                                                      stride_lightfield);
}

Transform::~Transform()
{
    delete ctx;
}

void Transform::dct_4d(const float *input,
                       float *output,
                       const Point4D &size,
                       const Point4D &origSize)
{
    auto *size_arr = size.to_array();
    ctx->forward(input, output, size_arr);
    delete[] size_arr;
}

void Transform::idct_4d(const float *input,
                        float *output,
                        const Point4D &size,
                        const Point4D &origSize)
{
    auto *size_arr = size.to_array();
    ctx->inverse(input, output, size_arr);
    delete[] size_arr;
}
