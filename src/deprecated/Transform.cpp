
#include <cassert>
#include "Transform.h"

namespace old {

Transform::Transform(const Point4D &dimBlock) {
    this->alocate_coeff_dct(dimBlock);
    this->input_4D_t = new float[dimBlock.getNSamples()];
}

void Transform::foward(const LFSample *input, float *output) {

}

void Transform::inverse(const float *input, LFSample *output) {

}

float *Transform::generate_dct_coeff(int N) {
    auto *output = new float[N * N];
    auto *p_output = output;
    for (int k = 0; k < N; k++) {
        // Normalization of parameters
        double s = (k == 0) ? 1 / (double) sqrt(N) : (sqrt(((double) 2 / N)));

        for (int n = 0; n < N; n++)
            *p_output++ = s * (double) cos((double) M_PI * (2 * n + 1) * k / (2 * N));
    }

    return output;
}

void Transform::dct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize) {
    if (this->dim_block != size) {
        this->delete_coeff_dct();
        this->alocate_coeff_dct(size);
    }

    const uint size_x = size.x,
            size_y = size.y,
            size_u = size.u,
            size_v = size.v;

    // TODO: Refactor these variables to something more meaninful 
    const uint offset1 = origSize.x,
            offset2 = offset1 * origSize.y,
            offset3 = offset2 * origSize.u,
            offset4 = offset3 * origSize.v;

    // input_4D_t[:] = input[:]
    std::copy(input, input + offset4, input_4D_t);
    float *p_input = input_4D_t, *p_output = output;

    for (int v = 0; v < size_v; ++v) { // TODO: positions
        for (int u = 0; u < size_u; ++u) {
            for (int y = 0; y < size_y; ++y) {
                dct_1D(p_input, p_output, this->coeff_dct1D[0], 1, size_x); // X AXIS

                p_input += offset1;
                p_output += offset1;
            }
            p_input += (origSize.y - size.y) * offset1;
            p_output += (origSize.y - size.y) * offset1;
        }
        p_input += (origSize.u - size.u) * offset2;
        p_output += (origSize.u - size.u) * offset2;
    }

    // input_4D_t[:] = output[:]
    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int v = 0; v < size_v; ++v) {
        for (int u = 0; u < size_u; ++u) {
            for (int x = 0; x < size_x; ++x) {
                dct_1D(p_input, p_output, this->coeff_dct1D[1], offset1, size_y); // Y AXIS

                ++p_input, ++p_output;
            }
            p_input += (offset2 - offset1) + (origSize.x - size.x);
            p_output += (offset2 - offset1) + (origSize.x - size.x);
        }
        p_input += (origSize.u - size.u) * offset2;
        p_output += (origSize.u - size.u) * offset2;
    }

    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int v = 0; v < size_v; ++v) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                dct_1D(p_input, p_output, this->coeff_dct1D[2], offset2, size_u); // U AXIS

                ++p_input, ++p_output;
            }
            p_input += (origSize.x - size.x);
            p_output += (origSize.x - size.x);
        }
        p_input += (offset3 - offset2) + (origSize.y - size.y) * offset1;
        p_output += (offset3 - offset2) + (origSize.y - size.y) * offset1;
    }

    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int u = 0; u < size_u; ++u) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                dct_1D(p_input, p_output, this->coeff_dct1D[3], offset3, size_v); // V AXIS

                ++p_input, ++p_output;
            }
            p_input += (origSize.x - size.x);
            p_output += (origSize.x - size.x);
        }
        p_input += (origSize.y - size.y) * offset1;
        p_output += (origSize.y - size.y) * offset1;
    }
}



/*********** Algorithm 
 * void dct(float[] in, float[] out, float[][] coeff, const size) {
 *      double sum;
 *      for (int k = 0; k < size; ++k) {
 *          sum = 0;
 *          for (int n = 0; n < size; ++n)
 *              sum += in[n] * coeff[k][n];
 *          out[k] = sum;
 *      }
 * }
 */
void Transform::dct_1D(const float *in, float *out, float *coeff, const uint offset, const uint size) {
    int N = size;

    float *p_table = coeff, *p_in, *p_out = out;
    double sum;

    for (int k = 0; k < N; ++k) {
        sum = 0;
        p_in = (float *) in;

        for (int n = 0; n < N; ++n) {
            sum += *p_in * *p_table;

            ++p_table;
            p_in += offset;
        }

        *p_out = (float) sum;
        p_out += offset;
    }

}

/*********** Algorithm 
 * void idct(float[] in, float[] out, float[][] coeff, const size) {
 *      double sum;
 *      for (int k = 0; k < size; ++k) {
 *          sum = 0;
 *          for (int n = 0; n < size; ++n)
 *              sum += in[n] * coeff[n][k];
 *          out[k] = sum;
 *      }
 * }
 */
void Transform::idct_1D(const float *in, float *out, float *coeff, const uint offset, const uint size) {
    int N = size;

    float *p_table, *p_in, *p_out = out;
    float sum;
    float partial_sum[size];

    for (int k = 0; k < N; ++k) {
        sum = 0;
        p_in = (float *) in;
        p_table = &coeff[k];
        for (int n = 0; n < N; ++n) {
            sum += *p_in * *p_table;
            partial_sum[n] = sum;

            p_table += N;
            p_in += offset;
        }

        *p_out = sum;
        p_out += offset;
    }

}

void Transform::idct_4d(const float *input, float *output, const Point4D &size, const Point4D &origSize) {
    if (this->dim_block != size) {
        this->delete_coeff_dct();
        this->alocate_coeff_dct(size);
    }

    const uint size_x = size.x,
            size_y = size.y,
            size_u = size.u,
            size_v = size.v;

    const uint offset1 = origSize.x,
            offset2 = offset1 * origSize.y,
            offset3 = offset2 * origSize.u,
            offset4 = offset3 * origSize.v;

    std::copy(input, input + offset4, input_4D_t);
    float *p_input = input_4D_t, *p_output = output;

    for (int u = 0; u < size_u; ++u) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                idct_1D(p_input, p_output, this->coeff_dct1D[3], offset3, size_v); // V AXIS

                ++p_input, ++p_output;
            }
            p_input += (origSize.x - size.x);
            p_output += (origSize.x - size.x);
        }
        p_input += (origSize.y - size.y) * offset1;
        p_output += (origSize.y - size.y) * offset1;
    }

    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int v = 0; v < size_v; ++v) {
        for (int y = 0; y < size_y; ++y) {
            for (int x = 0; x < size_x; ++x) {
                idct_1D(p_input, p_output, this->coeff_dct1D[2], offset2, size_u); // U AXIS

                ++p_input, ++p_output;
            }
            p_input += (origSize.x - size.x);
            p_output += (origSize.x - size.x);
        }
        p_input += (offset3 - offset2) + (origSize.y - size.y) * offset1;
        p_output += (offset3 - offset2) + (origSize.y - size.y) * offset1;
    }

    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int v = 0; v < size_v; ++v) {
        for (int u = 0; u < size_u; ++u) {
            for (int x = 0; x < size_x; ++x) {
                idct_1D(p_input, p_output, this->coeff_dct1D[1], offset1, size_y); // Y AXIS

                ++p_input, ++p_output;
            }
            p_input += (offset2 - offset1) + (origSize.x - size.x);
            p_output += (offset2 - offset1) + (origSize.x - size.x);
        }
        p_input += (origSize.u - size.u) * offset2;
        p_output += (origSize.u - size.u) * offset2;
    }

    std::copy(output, output + offset4, input_4D_t);
    p_input = input_4D_t, p_output = output;

    for (int v = 0; v < size_v; ++v) {
        for (int u = 0; u < size_u; ++u) {
            for (int y = 0; y < size_y; ++y) {
                idct_1D(p_input, p_output, this->coeff_dct1D[0], 1, size_x); // X AXIS

                p_input += offset1;
                p_output += offset1;
            }
            p_input += (origSize.y - size.y) * offset1;
            p_output += (origSize.y - size.y) * offset1;
        }
        p_input += (origSize.u - size.u) * offset2;
        p_output += (origSize.u - size.u) * offset2;
    }
}

void Transform::delete_coeff_dct() {
    delete[] this->coeff_dct1D[0];

    if (this->dim_block.y != this->dim_block.x) {
        delete[] this->coeff_dct1D[1];

        if (this->dim_block.u != this->dim_block.y) {
            delete[] this->coeff_dct1D[2];

            if (this->dim_block.v != this->dim_block.u) {
                delete[] this->coeff_dct1D[3];
            }
        }
    }
}

Transform::~Transform() {
    delete[] this->input_4D_t;

    this->delete_coeff_dct();
}

void Transform::alocate_coeff_dct(const Point4D &dimBlock) {
    this->dim_block = dimBlock;
    this->coeff_dct1D[0] = generate_dct_coeff(this->dim_block.x);

    this->coeff_dct1D[1] =
            (this->dim_block.y == this->dim_block.x) ? this->coeff_dct1D[0] : generate_dct_coeff(this->dim_block.y);

    this->coeff_dct1D[2] =
            (this->dim_block.u == this->dim_block.x) ? this->coeff_dct1D[0] :
            (this->dim_block.u == this->dim_block.y) ? this->coeff_dct1D[1] : generate_dct_coeff(this->dim_block.u);

    this->coeff_dct1D[3] =
            (this->dim_block.v == this->dim_block.x) ? this->coeff_dct1D[0] :
            (this->dim_block.v == this->dim_block.y) ? this->coeff_dct1D[1] :
            (this->dim_block.v == this->dim_block.u) ? this->coeff_dct1D[2] : generate_dct_coeff(this->dim_block.v);
}
};