
#ifndef QUANTIZATION_H
#define QUANTIZATION_H

#include <cmath>
#include "Typedef.h"
#include "Point4D.h"

class Quantization {
private:
    Point4D dim_block, weight_100;

    float *volume{nullptr};

    float qp;

    static float clip_min(float val, float min);

    void generateVolume();

public:
    Quantization(const Point4D &dimBlock, float qp, const Point4D &weight_100);

    void delete_volume();

    virtual ~Quantization();

    void inverse(const float *input, float *output);

    void foward(const float *input, float *output);
};


#endif //QUANTIZATION_H
