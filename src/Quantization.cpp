#include "Quantization.h"
#include "utils.h"
#include <cmath>

std::unordered_map<const Point4D, float *, Point4DHasher> Quantization::cache_lee;
std::unordered_map<const Point4D, float *, Point4DHasher> Quantization::cache_hayan;

Quantization::Quantization(const Point4D &dimBlock, const EncoderParameters codec_parameters):
    dim_block(dimBlock) {
    qp = codec_parameters.getQp();
    weight_100 = codec_parameters.quant_weight_100;
    lee_c = codec_parameters.lee_c;
    lee_ai = codec_parameters.lee_ai;
    lee_a0 = codec_parameters.lee_a0;
    lee_bi = codec_parameters.lee_bi;
    lee_b0 = codec_parameters.lee_b0;
}

void Quantization::foward(const float *input, float *output) { foward(HAIYAN, input, output); }

void Quantization::inverse(const float *input, float *output) { inverse(HAIYAN, input, output); }

float *Quantization::get_volume(const QuantizationVolumeType type) {
    auto cache = &Quantization::cache_hayan;
    float *volume = nullptr;
    switch (type) {
        case LEE: cache = &cache_lee;
        default: break;
    }
    try {
        volume = cache->at(dim_block);
    } catch (...) {
        volume = generate_volume(type);
        (*cache)[dim_block] = volume;
    }
    return volume;
}
void Quantization::set_lee_parameters(float C, float ai, float a0, float bi, float b0)
{
    lee_c = C; 
    lee_ai = ai; 
    lee_a0 = a0; 
    lee_bi = bi;
    lee_b0 = b0;
}

void Quantization::inverse(const QuantizationVolumeType type, const float *input, float *output) {
    const uint i_lim = this->dim_block.getNSamples();
    float *volume = get_volume(type);

    for (uint i = 0; i < i_lim; ++i)
        output[i] = (((float)input[i]) * volume[i]);
}

void Quantization::foward(const QuantizationVolumeType type, const float *input, float *output) {
    const uint i_lim = this->dim_block.getNSamples();
    float *volume = get_volume(type);
    for (uint i = 0; i < i_lim; ++i)
        output[i] = std::trunc(input[i] / volume[i]);
}

float Quantization::clip_min(float val, float min) { return (val > min) ? val : min; }



float *Quantization::generate_volume(const QuantizationVolumeType type) {
    const float weight_x = weight_100.x / 100.0f;
    const float weight_y = weight_100.y / 100.0f;
    const float weight_u = weight_100.u / 100.0f;
    const float weight_v = weight_100.v / 100.0f;

    float *volume = new float[dim_block.getNSamples()];
    float val;
    auto *it_vol = volume;
    volatile float qx, qy, qu, qv;

    for (int v = 0; v < dim_block.v; ++v) {
        qv = v * weight_v;
        for (int u = 0; u < dim_block.u; ++u) {
            qu = u * weight_u;
            for (int y = 0; y < dim_block.y; ++y) {
                qv = y * weight_y;
                for (int x = 0; x < dim_block.x; ++x) {
                    qx = x * weight_x;
                    switch (type) {
                        case LEE: {
                            float prod = (qx + 1.0F) * (qy + 1.0F) * (qu + 1.0F) * (qv + 1.0F);
                            if (qx * qy * qu * qv <= lee_c)
                                val = 1 + lee_ai * (1 - std::exp(-lee_bi * prod) / std::exp(-lee_bi));
                            else
                                val = lee_a0 * (1 - std::exp(-lee_b0 * prod));
                        } break;
                        case HAIYAN: {
                            val = 4 + qx + qy + qu + qv;
                        } break;
                    }
                    val *= qp;
                    *it_vol++ = val > 1 ? val : 1.0F;
                }
            }
        }
    }
    if (type == HAIYAN) volume[0] = 1.0;
    return volume;
}

void Quantization::delete_volume() {}

Quantization::~Quantization() { this->delete_volume(); }
