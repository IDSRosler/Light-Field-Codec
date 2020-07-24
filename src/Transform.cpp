
#include "Transform.h"
#include "Quantization.h"
#include "EncBitstreamWriter.h"
#include "LRE.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include "utils.h"

std::map<size_t, float *> Transform::cache_dct_ii;
std::map<size_t, float *> Transform::cache_dst_i;
std::map<size_t, float *> Transform::cache_dst_vii;

float *Transform::get_coefficients(Transform::TransformType type, const size_t size) {
    float *coeff = nullptr;
    switch (type) {
        case DST_I: try { coeff = cache_dst_i.at(size);
            } catch (...) {
                coeff = sd_dst_i(size);
                cache_dst_i[size] = coeff;
            }
            break;
        case DST_VII: try { coeff = cache_dst_vii.at(size);
            } catch (...) {
                coeff = sd_dst_vii(size);
                cache_dst_vii[size] = coeff;
            }
            break;
        case DCT_II:
        default: try { coeff = cache_dct_ii.at(size);
            } catch (...) {
                coeff = sd_dct_ii(size);
                cache_dct_ii[size] = coeff;
            }
    }
    return coeff;
}

void Transform::sd_forward(Transform::TransformType type,
                           const float *in,
                           float *out,
                           const size_t offset,
                           const size_t size) {
    float *coeff = get_coefficients(type, size);
    auto pout = out;
    auto pcoeff = coeff;
    for (int k = 0; k < size; k++, pout += offset) {
        auto pin = in;
        *pout = 0;
        for (int n = 0; n < size; n++, pin += offset, pcoeff++)
            *pout += *pin * *pcoeff;
    }
}

void Transform::sd_inverse(Transform::TransformType type,
                           const float *in,
                           float *out,
                           const size_t offset,
                           const size_t size) {
    float *coeff = get_coefficients(type, size);
    auto pout = out;
    for (int k = 0; k < size; k++, pout += offset) {
        auto pin = in;
        auto pcoeff = coeff + k;
        *pout = 0;
        for (int n = 0; n < size; n++, pin += offset, pcoeff += size)
            *pout += *pin * *pcoeff;
    }
}

Transform::Transform(Point4D &shape) {
    this->shape = shape;
    this->enforce_transform = TransformType::NO_TRANSFORM;
    stride = make_stride(shape);
    flat_size = stride.v * shape.v;
    partial_values = new float[flat_size];
    fake_encoder = nullptr;    
}

Transform::Transform(EncoderParameters &params) : Transform(params.dim_block)
{
    codec_parameters = params;
    fake_encoder = new EncBitstreamWriter(&codec_parameters, 10'000'000, true);
}

Transform::~Transform()
{
    delete fake_encoder;
    delete[] partial_values;
}

float *Transform::sd_dst_vii(size_t size) {
    auto *output = new float[size * size];
    auto *pout = output;
    for (int k = 0; k < size; k++) {
        double s = (double)2 / (sqrt(2 * size + 1));
        for (int n = 0; n < size; n++)
            *pout++ = s * sin((double)M_PI * (n + 1) * (2 * k + 1) / (2 * size + 1));
    }
    return output;
}
float *Transform::sd_dst_i(size_t size) {
    auto *output = new float[size * size];
    auto *pout = output;
    double s = sqrt(2.0L / (size + 1.0L));
    for (int k = 0; k < size; k++) {
        for (int n = 0; n < size; n++) {
            double theta = (M_PI * (k + 1.0L) * (n + 1.0L)) / (size + 1.0L);
            *pout++ = s * sin(theta);
        }
    }
    return output;
}

float *Transform::sd_dct_ii(size_t size) {
    float *output = new float[size * size];
    float *pout = output;
    for (int k = 0; k < size; k++) {
        double s = (k == 0) ? 1 / (double)sqrt(size) : (sqrt(((double)2 / size)));
        for (int n = 0; n < size; n++)
            *pout++ = s * cos((double)M_PI * (2 * n + 1) * k / (2 * size));
    }
    return output;
}

auto Transform::get_transform_vector(TransformType transform) {
    switch (transform) {
        case HYBRID: return std::vector({DCT_II, DCT_II, DST_I, DST_I});
        default: return std::vector({transform, transform, transform, transform});
    }
}


void Transform::md_fw_axis(int ax, TransformType type, float *input, float *output, Point4D &shape) {
    int axis[3] = {0};
    auto *p = axis;
    for (int x = 0; x < 4; x++)
        if (x != ax)
            *p++ = x;
    volatile int index;
    for (int k = 0; k < shape[axis[2]]; ++k)
        for (int j = 0; j < shape[axis[1]]; ++j)
            for (int i = 0; i < shape[axis[0]]; ++i) {
                index = k * stride[axis[2]] +
                        j * stride[axis[1]] + 
                        i * stride[axis[0]];
                sd_forward(type,
                           input + index, 
                           output + index, 
                           stride[ax], 
                           shape[ax]);
            }
}

void Transform::md_in_axis(int ax, TransformType type, float *input, float *output, Point4D &shape) {
    int axis[3] = {0};
    auto *p = axis;
    for (int x = 0; x < 4; x++)
        if (x != ax)
            *p++ = x;
    volatile int index;
    for (int k = 0; k < shape[axis[2]]; ++k)
        for (int j = 0; j < shape[axis[1]]; ++j)
            for (int i = 0; i < shape[axis[0]]; ++i) {
                index = k * stride[axis[2]] +
                        j * stride[axis[1]] + 
                        i * stride[axis[0]];
                sd_inverse(type,
                           input + index, 
                           output + index, 
                           stride[ax], 
                           shape[ax]);
            }
}

void Transform::md_forward(TransformType type, float *input, float *output, Point4D &shape)
{
    int ax_order[] = {0, 1, 2, 3};
    float *pin = input;
    float *pout = output;

    for (int i = 0; i < 4; i++) {
        md_fw_axis(ax_order[i], type, pin, pout, shape);
        pin = partial_values, pout = output;
        if (i < 3)
            std::copy(output, output + flat_size, partial_values);
    }
}

void Transform::md_inverse(TransformType type, float *input, float *output, Point4D &shape) {
    int ax_order[] = {3, 2, 1, 0};

    float *pin = (float *)input;
    float *pout = output;

    for (int i = 0; i < 4; i++) {
        md_in_axis(ax_order[i], type, pin, pout, shape);
        pin = partial_values, pout = output;
        if (i < 3)
            std::copy(output, output + flat_size, partial_values);
    }
}

void Transform::flush_cache() {}

Transform::TransformType Transform::get_type(std::string transform) {
    if (transform == "DST")
        return DST_I;
    else if (transform == "DST_I")
        return DST_I;
    else if (transform == "DST_VII")
        return DST_VII;
    else if (transform == "DCT")
        return DCT_II;
    else if (transform == "MULTI")
        return MULTI;
    else if (transform == "HYBRID")
        return HYBRID;
    else
        std::cerr << "Unkown transform: " << transform << std::endl;
    return DCT_II;
}



void Transform::set_position(int channel, const Point4D &current_pos) {
    this->channel = channel;
    this->position = current_pos;
}

auto Transform::get_quantization_procotol(TransformType transform)
{
    return Quantization::HAIYAN;
    switch (transform)
    {
    case DCT_II:
        return Quantization::LEE;
    default:
        return Quantization::HAIYAN;
    }
}
static float mse(float *original, float *reconstructed, unsigned size) {
    auto value = 0;
    for (int i = 0; i < size; i++)
        value += std::pow(original[i] - reconstructed[i], 2);
    return value / size;
}

auto Transform::calculate_distortion(float *block, float *result, Point4D &shape) {
    shape.updateNSamples();
    const int SIZE = this->shape.getNSamples();
    std::vector<TransformType> transforms;

    if (enforce_transform == NO_TRANSFORM)
        return std::tuple(NO_TRANSFORM, std::numeric_limits<float>::infinity());
    // TODO: Enable transforms by macro
    else if (enforce_transform == MULTI)
        for (const auto &type: ALL_TRANSFORMS)
            transforms.push_back(type);
    else
        transforms.push_back(enforce_transform);

    float blk_transf[SIZE];
    float blk_qnt[SIZE];
    float blk_iquant[SIZE];
    float blk_itransf[SIZE];

    std::fill(blk_transf, blk_transf + SIZE, 0);
    std::fill(blk_qnt, blk_qnt + SIZE, 0);
    std::fill(blk_iquant, blk_iquant + SIZE, 0);
    std::fill(blk_itransf, blk_itransf + SIZE, 0);

    float best_score = std::numeric_limits<float>::infinity();
    TransformType best_type;

    Quantization quantizer(shape, codec_parameters);
    for (const auto &type: transforms) {
        md_forward(type, block, blk_transf, shape);
        quantizer.foward(get_quantization_procotol(type), blk_transf, blk_qnt);
        quantizer.inverse(get_quantization_procotol(type), blk_qnt, blk_iquant);
        md_inverse(type, blk_iquant, blk_itransf, shape);
        auto curr_score = mse(block, blk_itransf, SIZE);

        if (curr_score < best_score) {
            best_type = type;
            best_score = curr_score;
            std::copy(blk_qnt, blk_qnt + SIZE, result);
        }
    }
    return std::tuple(best_type, best_score);
}

auto Transform::calculate_tree(float *block, float *result, Point4D &shape, int level) {
    if (!block) /* Dummy return to deduce return type */
        return std::tuple((Node *)nullptr, 0.0F);

    const int SIZE = this->shape.getNSamples();
    shape.updateNSamples();

    float block_transformed[SIZE];
    float block_segmented[SIZE];
    float segments_results[SIZE];
    float joined_results[SIZE];

    std::fill(block_transformed, block_transformed + SIZE, 0);
    std::fill(block_segmented, block_segmented + SIZE, 0);
    std::fill(segments_results, segments_results + SIZE, 0);
    std::fill(joined_results, joined_results + SIZE, 0);

    float segments_score = std::numeric_limits<float>::infinity();
    auto block_shape = this->shape.to_vector();
    float *transformed_block;
    float final_score;

    auto [block_type, block_score] = calculate_distortion(block, block_transformed, shape);
    auto root = new Node();

    if (level > 0) {
        float segments_score_sum = 0;
        auto shapes = make_shapes(block_shape, 1);

        auto shape_bak = this->shape;
        auto stride_bak = this->stride;
        auto flat_size_bak = this->flat_size;
        segment_block(block, block_shape, block_segmented, 1);
        float *pseg = block_segmented;
        float *pres = segments_results;

        for (int i = 0; i < shapes.size(); i++) {
            Point4D seg_shape(shapes[i].data());
            this->shape = seg_shape;
            this->flat_size = seg_shape.getNSamples();
            this->stride = make_stride(seg_shape);
            auto [node, score] = calculate_tree(pseg, pres, seg_shape, level - 1);

            root->set_child(i, node);
            segments_score_sum += score;
            pseg += seg_shape.getNSamples();
            pres += seg_shape.getNSamples();
        }
        this->shape = shape_bak;
        this->stride = stride_bak;
        this->flat_size = flat_size_bak;

        segments_score = segments_score_sum / shapes.size();
    }
    if (segments_score < block_score) {
        transformed_block = segments_results;
        final_score = segments_score;
    } else {
        root->set_transform_type(block_type);
        transformed_block = block_transformed;
        final_score = block_score;
    }
    std::copy(transformed_block, transformed_block + SIZE, result);
    return std::tuple(root, final_score);
}

void Transform::reconstruct_from_tree(Node *root, float *input, float *output, Point4D &shape) {
    shape.updateNSamples();
    const int SIZE = this->shape.getNSamples();
    float block_iquant[SIZE];
    float block_segments[SIZE];
    float segment_results[SIZE];
    std::fill(block_iquant, block_iquant + SIZE, 0);
    std::fill(block_segments, block_segments + SIZE, 0);
    std::fill(segment_results, segment_results + SIZE, 0);

    if (root->transform_type != -1) {
        Quantization quantizer(shape, codec_parameters);
        quantizer.inverse(get_quantization_procotol((TransformType)root->transform_type), input,
                          block_iquant);
        md_inverse((TransformType)root->transform_type, block_iquant, output, shape);
    } else {
        auto block_shape = this->shape.to_vector();
        auto shapes = make_shapes(block_shape, 1);
        auto shape_bak = this->shape;
        auto stride_bak = this->stride;
        auto flat_size_bak = this->flat_size;

        float *pseg = input;
        float *pres = segment_results;
        for (int i = 0; i < shapes.size(); i++) {
            Point4D seg_shape(shapes[i].data());
            this->shape = seg_shape;
            this->flat_size = seg_shape.getNSamples();
            this->stride = make_stride(seg_shape);
            reconstruct_from_tree(root->children[i], pseg, pres, seg_shape);
            pseg += seg_shape.getNSamples();
            pres += seg_shape.getNSamples();
        }
        this->shape = shape_bak;
        this->stride = stride_bak;
        this->flat_size = flat_size_bak;
        join_segments(segment_results, block_shape, output, 1);
    }
}

std::string
Transform::forward(TransformType transform, float *input, float *output, Point4D &shape) {
    shape.updateNSamples();
    char buffer[1 << 16];
#if !!LFCODEC_FORCE_DCT_NON_LUMA && USE_YCbCr == 1
    this->enforce_transform = channel == 0 ? transform : DCT_II;
#else
    this->enforce_transform = transform;
#endif
    const int MAX_LEVELS = disable_segmentation ? 0 : std::min(LFCODEC_SEGMENTATION_MAX_LEVELS, 3);
#if !!LFCODEC_USE_SEGMENTATION
#if !!LFCODEC_FORCE_DCT_NON_LUMA && USE_YCbCr == 1
    std::vector channel_mapping = {MAX_LEVELS, 0, 0};
#else
    std::vector channel_mapping = {MAX_LEVELS, MAX_LEVELS, MAX_LEVELS};
#endif
#else
    std::vector channel_mapping = {0, 0, 0};
#endif
    auto [root, mse] = calculate_tree(input, output, shape, channel_mapping[channel]);
    root->to_string(buffer);
    std::string tree = buffer;
    delete root;
    return tree;
}

void Transform::inverse(TransformType transform, float *input, float *output, Point4D &shape) {
    shape.updateNSamples();
    const int SIZE = this->shape.getNSamples();
    float block_iquant[SIZE];
    Quantization quantizer(shape, codec_parameters);
    quantizer.inverse(get_quantization_procotol(transform), input, block_iquant);
    md_inverse(transform, block_iquant, output, shape);
}

void Transform::inverse(const std::string tree, float *input, float *output, Point4D &shape) {
    Node root;
    root.from_string(tree.c_str());
    reconstruct_from_tree(&root, input, output, shape);
}
