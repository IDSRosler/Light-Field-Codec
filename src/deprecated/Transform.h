
#ifndef _OLD_TRANSFORM_H
#define _OLD_TRANSFORM_H

#include <zconf.h>
#include <cmath>
#include "Typedef.h"
#include "../Point4D.h"

namespace old {

class Transform {
private:
    float *coeff_dct1D[4]{};
    Point4D dim_block;

    float *input_4D_t;

    static float *generate_dct_coeff(int N);

    static void dct_1D(const float *in, float *out, float *coeff, const uint offset, const uint size);

    static void idct_1D(const float *in, float *out, float *coeff, const uint offset, const uint size);


public:
    explicit Transform(const Point4D &dimBlock);

    void dct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize);

    void idct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize);

    ~Transform();

    void delete_coeff_dct();

    void alocate_coeff_dct(const Point4D &dimBlock);

    void foward(const short *input, float *output);

    void inverse(const float *input, short *output);
};
};

#endif //TRANSFORM_H
