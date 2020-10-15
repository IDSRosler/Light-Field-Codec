
#include <string>
#include "Transform.h"
#include "Quantization.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>


Transform::Transform(Point4D &shape) {
    this->block_shape = shape;
    block_stride = make_stride(shape);
    flat_size = block_stride.v * shape.v;
    m_partial_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_r_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_tf_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_lre_block = std::unique_ptr<int[]>(new int[flat_size]);
    lre = std::make_unique<LRE>(shape);
}

Transform::Transform(EncoderParameters &params) : Transform(params.dim_block) {
    codec_parameters = params;
    fake_encoder = std::make_unique<EncBitstreamWriter>(&codec_parameters, 1048576, true);
}

float *Transform::get_coefficients(Transform::TransformType type,
                                   const size_t size) {
    float *coeff;
    switch (type) {
        case DST_I:
            try {
                coeff = cache_dst_i.at(size);
            }
            catch (...) {
                coeff = sd_dst_i(size);
                cache_dst_i[size] = coeff;
            }
            break;
        case DST_VII:
            try {
                coeff = cache_dst_vii.at(size);
            }
            catch (...) {
                coeff = sd_dst_vii(size);
                cache_dst_vii[size] = coeff;
            }
            break;
        case DCT_II:
        default:
            try {
                coeff = cache_dct_ii.at(size);
            }
            catch (...) {
                coeff = sd_dct_ii(size);
                cache_dct_ii[size] = coeff;
            }
    }
    return coeff;
}


void Transform::sd_forward(Transform::TransformType type, const float *in,
                           float *out, const size_t stride, const size_t size) {
    auto naive_approach = [&]() {
        float *coeff = get_coefficients(type, size);
        auto pcoeff = coeff;
        auto pout = out;
        for (std::size_t k = 0; k < size; k++, pout += stride) {
            auto pin = in;
            *pout = 0;
            for (std::size_t n = 0; n < size; n++, pin += stride, pcoeff++)
                *pout += *pin * *pcoeff;
        }
    };

    naive_approach();

}

void Transform::sd_inverse(Transform::TransformType type, const float *in,
                           float *out, const size_t stride, const size_t size) {
    auto naive_approach = [&]() {
        float *coeff = get_coefficients(type, size);
        auto pout = out;
        for (std::size_t k = 0; k < size; k++, pout += stride) {
            auto pin = in;
            auto pcoeff = coeff + k;
            *pout = 0;
            for (std::size_t n = 0; n < size; n++, pin += stride, pcoeff += size)
                *pout += *pin * *pcoeff;
        }
    };

    naive_approach();

}

float *Transform::sd_dst_vii(size_t size) {
    auto *output = new float[size * size];
    auto *pout = output;
    for (std::size_t k = 0; k < size; k++) {
        double s = (double) 2 / (sqrt(2 * size + 1));
        for (std::size_t n = 0; n < size; n++)
            *pout++ = s * sin((double) M_PI * (n + 1.0) * (2.0 * k + 1) / (2.0 * size + 1));
    }
    return output;
}

float *Transform::sd_dst_i(size_t size) {
    auto *output = new float[size * size];
    auto *pout = output;
    double s = sqrt(2.0 / (size + 1.0));
    for (std::size_t k = 0; k < size; k++) {
        for (std::size_t n = 0; n < size; n++) {
            double theta = (M_PI * (k + 1.0) * (n + 1.0)) / (size + 1.0);
            *pout++ = s * sin(theta);
        }
    }
    return output;
}

float *Transform::sd_dct_ii(size_t size) {
    auto *output = new float[size * size];
    float *pout = output;
    for (std::size_t k = 0; k < size; k++) {
        double s = (k == 0) ? (1.0 / sqrt(size)) : (sqrt((2.0 / size)));
        for (std::size_t n = 0; n < size; n++)
            *pout++ = static_cast<float>(s * cos((double) M_PI * ((2.0 * n + 1) * k) / (2.0 * size)));
    }
    return output;
}


void Transform::md_forward_single_axis(const int ax,
                                       const TransformType type,
                                       const float *input,
                                       float *output,
                                       const Point4D &shape) {
    using index_type = Point4D::value_type;

    index_type ax_shape = shape[ax];
    index_type ax_stride = block_stride[ax];

    int axis[3] = {0};
    auto *p = axis;


    for (int x = 0; x < 4; x++)
        if (x != ax)
            *p++ = x;

    std::vector<index_type> i_indexes(shape[axis[0]]);
    std::vector<index_type> j_indexes(shape[axis[1]]);
    std::vector<index_type> k_indexes(shape[axis[2]]);

    std::iota(std::begin(i_indexes), std::end(i_indexes), 0);
    std::iota(std::begin(j_indexes), std::end(j_indexes), 0);
    std::iota(std::begin(k_indexes), std::end(k_indexes), 0);

    std::for_each(std::cbegin(k_indexes), std::cend(k_indexes), [&](index_type k) {
        std::for_each(std::cbegin(j_indexes), std::cend(j_indexes), [&](index_type j) {
            std::for_each(std::cbegin(i_indexes), std::cend(i_indexes), [&](index_type i) {
                index_type index = k * block_stride[axis[2]] +
                                   j * block_stride[axis[1]] +
                                   i * block_stride[axis[0]];
                sd_forward(type, input + index, output + index, ax_stride, ax_shape);
            });
        });
    });


}

void Transform::md_inverse_single_axis(const int ax,
                                       const TransformType type,
                                       const float *input,
                                       float *output,
                                       const Point4D &shape) {

    using index_type = Point4D::value_type;
    int axis[3] = {0};
    auto *p = axis;
    for (int x = 0; x < 4; x++)
        if (x != ax)
            *p++ = x;

    std::vector<index_type> i_indexes(shape[axis[0]]);
    std::vector<index_type> j_indexes(shape[axis[1]]);
    std::vector<index_type> k_indexes(shape[axis[2]]);
    std::iota(std::begin(i_indexes), std::end(i_indexes), 0);
    std::iota(std::begin(j_indexes), std::end(j_indexes), 0);
    std::iota(std::begin(k_indexes), std::end(k_indexes), 0);
    index_type ax_shape = shape[ax];
    index_type ax_stride = block_stride[ax];

    std::for_each(std::cbegin(k_indexes), std::cend(k_indexes), [&](index_type k) {
        std::for_each(std::cbegin(j_indexes), std::cend(j_indexes), [&](index_type j) {
            std::for_each(std::cbegin(i_indexes), std::cend(i_indexes), [&](index_type i) {
                index_type index = k * block_stride[axis[2]] +
                                   j * block_stride[axis[1]] +
                                   i * block_stride[axis[0]];
                sd_inverse(type, input + index, output + index, ax_stride, ax_shape);
            });
        });
    });
}


void Transform::md_forward(const TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &offset_,
                           const Point4D &shape) {
    int ax_order[] = {3, 2, 1, 0};
    float *partial_values = m_partial_block.get();
    Point4D adjusted_shape = shape;

    for (int ax = 0; ax < 4; ax++)
        if (offset_[ax] + shape[ax] > block_shape[ax])
            adjusted_shape[ax] = block_shape[ax] - offset_[ax];


    int offset = calc_offset(offset_, block_stride);
    const float *pin = input + offset;
    float *pout = output + offset;

    for (auto i : {0, 1, 2, 3}) {
        md_forward_single_axis(ax_order[i], type, pin, pout, adjusted_shape);
        pin = partial_values + offset;
        pout = output + offset;
        if (i < 3)
            // TODO: Optimize for smaller sizes
            std::copy(output, output + flat_size, partial_values);
    }


    if (codec_parameters.initialized && !codec_parameters.lossless) {
        if (codec_parameters.uniform_quantization) {
            FOREACH_4D_IDX(i, shape, block_stride)
                pout[i] = std::trunc(pout[i] / codec_parameters.qp);
        } else {
            auto volume_stride = make_stride(block_shape);
            Quantization q(block_shape, codec_parameters);
            float *volume = q.get_volume(Quantization::HAIYAN);
            FOREACH_4D_IDX(i, shape, block_stride) {
                auto volume_offset = calc_offset(x, y, u, v, volume_stride);
                pout[i] = std::trunc(pout[i] / volume[volume_offset]);
            }
        }
    }
}

void Transform::md_inverse(const TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &offset_,
                           const Point4D &shape) {
    Point4D adjusted_shape = shape;
    int ax_order[] = {0, 1, 2, 3};

    float *intermediary_values = m_partial_block.get();


    for (int ax = 0; ax < 4; ax++) {
        if (offset_[ax] + shape[ax] > block_shape[ax])
            adjusted_shape[ax] = block_shape[ax] - offset_[ax];
    }

    int offset = calc_offset(offset_, block_stride);
    float *pout = output + offset;
    const float *pin;

    float block[flat_size];

    if (codec_parameters.initialized && !codec_parameters.lossless) {
        pin = block + offset;
        if (codec_parameters.uniform_quantization) {
            FOREACH_4D_IDX(i, shape, block_stride)
                block[i + offset] = input[i + offset] * codec_parameters.qp;
        } else {
            auto volume_stride = make_stride(block_shape);
            Quantization q(block_shape, codec_parameters);
            float *volume = q.get_volume(Quantization::HAIYAN);

            FOREACH_4D_IDX(i, shape, block_stride) {
                auto volume_offset = calc_offset(x, y, u, v, volume_stride);
                block[i + offset] = input[i + offset] * volume[volume_offset];
            }
        }
    }


    for (int i = 0; i < 4; i++) {
        md_inverse_single_axis(ax_order[i], type, pin, pout, adjusted_shape);
        pin = intermediary_values + offset;
        pout = output + offset;
        if (i < 3)
            // TODO: Optimize for smaller sizes
            std::copy(output, output + flat_size, intermediary_values);
    }
}

void Transform::flush_cache() {}


Transform::TransformType Transform::get_type(const std::string &transform) {
    if (transform == "DST" || transform == "DST_I")
        return DST_I;
    else if (transform == "DST_VII")
        return DST_VII;
    else if (transform == "DCT" || transform == "DCT_II")
        return DCT_II;
    else if (transform == "MULTI")
        return MULTI;
    else if (transform == "HYBRID")
        return HYBRID;
    else
        std::cerr << "Unknown transform: " << transform << std::endl;
    return DCT_II;
}

void Transform::set_position(int _channel, const Point4D &current_pos) {
    this->channel = _channel;
    this->position = current_pos;
}

void Transform::forward(const std::string &descriptor, const float *block, float *result, const Point4D &_shape) {
    for (const auto &desc: parse_descriptor(descriptor, _shape, true)) {
        const auto &[offset, shape, type] = desc;
        md_forward(static_cast<TransformType>(type), block, result, offset, shape);
    }
}

void Transform::inverse(const std::string &descriptor, const float *block, float *result, const Point4D &_shape) {
    for (const auto &desc: parse_descriptor(descriptor, _shape, true)) {
        const auto &[offset, shape, type] = desc;
        md_inverse(static_cast<TransformType>(type), block, result, offset, shape);
    }
}

auto Transform::find_min_rd_cost(const float *block,
                                 float *transformed_block,
                                 const Point4D &shape,
                                 const Point4D &offset) {
    auto rec_block = m_temp_r_block.get();
    auto lre_block = m_temp_lre_block.get();
    const std::size_t size = shape.getNSamples();
    quadtree::node_type node = quadtree::make_node();
    double min_cost = std::numeric_limits<double>::infinity();
    auto linear_offset = calc_offset(offset, block_stride);

    for (auto &transform : Transform::T_CHOICES) {
        auto lre_block_ptr = lre_block;
        md_forward(transform, block, transformed_block, offset, shape);
        md_inverse(transform, transformed_block, rec_block, offset, shape);
        double sse = 0;
        FOREACH_4D_IDX(i, shape, block_stride) {
            double error = block[i+linear_offset] - rec_block[i+linear_offset];
            sse += error * error;
            *lre_block_ptr++ = static_cast<int>(std::round(transformed_block[i]));
        }

        std::size_t encoding_size;
        auto czi_result = lre->encodeCZI(lre_block, 0, size);
        double mse = sse / size;
        if (fake_encoder != nullptr && !codec_parameters.fast) {
            encoding_size = fake_encoder->write4DBlock(lre_block, size, czi_result);
            fake_encoder->reset();
        } else {
            encoding_size = czi_result.size();
        }

        double rd_cost = mse + codec_parameters.lambda * encoding_size / static_cast<double>(size);
        if (rd_cost < min_cost) {
            min_cost = rd_cost;
            node->transform = transform;
        }
    }
    return std::make_tuple(min_cost, node);
}


auto Transform::min_partition_tree(const float *block,
                                   float *transformed_block,
                                   const Point4D &from_shape,
                                   const Point4D &from_offset,
                                   size_t partition_level) {
    auto[min_cost, node] = find_min_rd_cost(block, transformed_block, from_shape, from_offset);

    if (partition_level == 0)
        return std::make_tuple(min_cost, node);

    double spatial_cost = std::numeric_limits<double>::infinity();
    quadtree::node_type spatial_tree = quadtree::make_node();
    spatial_tree->partition = 's';

    double angular_cost = std::numeric_limits<double>::infinity();
    quadtree::node_type angular_tree = quadtree::make_node();
    angular_tree->partition = 'a';

    auto min_spatial = codec_parameters.transform_min_spatial_size;
    auto min_angular = codec_parameters.transform_min_angular_size;

    if (from_shape.x > min_spatial && from_shape.y > min_spatial) {
        spatial_cost = 0;
        for (size_t c = 0; c < 4U; c++) {
            const auto &[shape, offset] = partition_at(from_shape, from_offset, c, 's');
            const auto &[cost, node] = min_partition_tree(block, transformed_block, shape, offset,partition_level - 1);
            spatial_cost += cost;
            spatial_tree->set_child(c, node);
        }
        spatial_cost /= 4;
    }

    if (from_shape.u > min_angular && from_shape.v > min_angular) {
        angular_cost = 0;
        for (size_t c = 0; c < 4U; c++) {
            const auto &[shape, offset] = partition_at(from_shape, from_offset, c, 'a');
            const auto &[cost, node] = min_partition_tree(block, transformed_block, shape, offset,
                                                          partition_level - 1);
            angular_cost += cost;
            angular_tree->set_child(c, node);
        }
        angular_cost /= 4;
    }
    if (spatial_cost < min_cost) {
        min_cost = spatial_cost;
        node = spatial_tree;
    }
    if (angular_cost < min_cost) {
        min_cost = angular_cost;
        node = angular_tree;
    }
    return std::make_tuple(min_cost, node);
}


std::pair<std::string, double> Transform::forward(const float *block, float *result, const Point4D &_shape) {
    auto[cost, node] = min_partition_tree(block, result, _shape, Point4D(0, 0, 0, 0), PARTITION_TREE_DEPTH);
    std::string descriptor = node->repr();
    forward(descriptor, block, result, _shape);
    return std::make_pair(descriptor, cost);
}
// EOF