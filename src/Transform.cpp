
#include "Transform.h"
#include "EncBitstreamWriter.h"
#include "LRE.h"
#include "Quantization.h"
#include "ScanOrder.h"
#include "utils.h"
#include "cstdlib"
#if LFCODEC_FAST_DCT
#include "FastImplementations/FastDctFFT.h"
#endif

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>

// #include <execution>
Transform::Transform(Point4D &shape) {
    this->block_shape = shape;
    this->enforce_transform = TransformType::NO_TRANSFORM;
    block_stride = make_stride(shape);
    flat_size = block_stride.v * shape.v;
    m_partial_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_r_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_tf_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_lre_block = std::unique_ptr<int[]>(new int[flat_size]);
    lre = std::make_unique<LRE>(shape);

    for (const auto &N : QUADTREE_NODES_COUNT) {
        for (const auto &tree : generate_full_binary_trees(N)) {
            tree_repr_vector.push_back(tree->repr());
        }
    }

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

#if LFCODEC_FAST_DCT

#else
    naive_approach();
#endif
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
#if LFCODEC_FAST_DCT

#else
    naive_approach();
#endif
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
                m_forward_count_n2 += ax_shape * ax_shape;
                auto log_ax_shape = std::log2l(ax_shape);
                assert(log_ax_shape < ax_shape);
                assert(ax_shape > 0);
                m_forward_count_nlogn += static_cast<std::uint64_t>(ax_shape * log_ax_shape);
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
                m_inverse_count_n2 += ax_shape * ax_shape;
                auto log_ax_shape = std::log2l(ax_shape);
                assert(log_ax_shape < ax_shape);
                assert(ax_shape > 0);
                m_inverse_count_nlogn += ax_shape * log_ax_shape;

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
    const float *pin = input + offset;

    float block[flat_size];
    pin = block + offset;

    if (codec_parameters.initialized && !codec_parameters.lossless) {
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

std::pair<std::string, double> Transform::forward(const float *block, float *result, const Point4D &_shape) {
    float *tf_block = m_temp_tf_block.get();
    m_min_rd_cost = std::numeric_limits<double>::infinity();
    for (const auto &tree_repr : tree_repr_vector) {

        std::deque<std::pair<Point4D, Point4D>> stack;
        stack.push_front(std::make_pair(_shape, Point4D(0, 0, 0, 0)));
        forward_fast(tree_repr, 0, block, tf_block, stack);
//#if true
//        std::cout << "[log] Pos(x=" << position.x << ","
//                               "y=" << position.y << ","
//                               "u=" << position.u << ","
//                               "v=" << position.v << ","
//                               "ch=" << channel << ") "
//                        "Descriptor(tree=" << tree_repr <<", "
//                                   "desc=" << m_min_descriptor << m_encoding_type << ", "
//                                   "rd_cost=" << m_min_rd_cost << ")\n";
//        std::cout.flush();
//#endif
    }

    forward(m_min_descriptor, block, result, _shape);
    return std::make_pair(m_min_descriptor, m_min_rd_cost);
}
/***/
void Transform::calculate_rd_cost(const float *block, const std::string &descriptor) {

    std::string min_descriptor;
    // flat_size changes during runtime, but here we need the full block size
    size_t size = block_shape.v * block_stride.v;
    int *lre_block = m_temp_lre_block.get();
    float *t_block = m_temp_tf_block.get();
    float *r_block = m_temp_r_block.get();
    // std::fill(lre_block, lre_block + size, 0);
    assert(size == flat_size);
    double sse = std::inner_product(block,
                                    block + size,
                                    r_block,
                                    0.0,
                                    [](auto acc, auto val) { return acc + val; },
                                    [](auto blk, auto rec) { auto error = blk - rec; return error * error; });

    std::transform(t_block, t_block + size, lre_block, std::truncf);
    std::size_t lre_size, czi_size, encoding_size;
    char encoding_type;
    auto mse = sse / size;
    auto lre_result = lre->encodeLRE(lre_block, size);
    auto czi_result = lre->encodeCZI(lre_block, 0, size);

    if (fake_encoder != nullptr) {
        lre_size = fake_encoder->write4DBlock(lre_block, size, lre_result);
        czi_size = fake_encoder->write4DBlock(lre_block, size, czi_result);
        fake_encoder->reset();
    } else {
        lre_size = lre_result.size();
        czi_size = czi_result.size();
    }

    if (lre_size < czi_size) {
        encoding_size = lre_size;
        encoding_type = 'L';
    } else {
        encoding_size = czi_size;
        encoding_type = 'Z';
    }

    double rd_cost = mse + codec_parameters.lambda * encoding_size;
#if false
    std::cout << "[log] Pos(x=" << position.x << ","
                           "y=" << position.y << ","
                           "u=" << position.u << ","
                           "v=" << position.v << ","
                           "ch=" << channel << ") "
                        "Descriptor(text=" << descriptor << encoding_type << ", "
                                    "mse=" << mse << ", "
                                    "lre_size=" << lre_size << ", "
                                    "czi_size=" << czi_size << ", "
                                    "rd_cost=" << rd_cost << ")\n";
    std::cout.flush();
#endif
    if (rd_cost < m_min_rd_cost) {
        m_min_rd_cost = rd_cost;
        m_min_descriptor = descriptor;
        m_encoding_type = encoding_type;
    }
}

void Transform::forward_fast(
        const std::string &descriptor,
        std::size_t index,
        const float *block,
        float *result,
        std::deque<std::pair<Point4D, Point4D>> stack) {
    TIMER_TIC(m_timer_substr);
    std::string prefix(descriptor.substr(0, index > 0 ? index : 0));
    std::string suffix(descriptor.substr(index + 1, descriptor.size() - 1));
    TIMER_TOC(m_timer_substr);

    const auto symbol = descriptor[index];
    auto[shape, offset] = stack.front();
    shape.updateNSamples();
    stack.pop_front();

    Quantization q(shape, codec_parameters);
    const float *q_volume = q.get_volume(Quantization::HAIYAN);

    switch (symbol) {
        case 'P': {
            for (const auto type : P_CHOICES) {
                std::stringstream ss;
                ss << prefix << type << suffix;
                TIMER_TIC(m_timer_split);
                auto temp_stack = split_coordinate(type, offset, shape);
                TIMER_TOC(m_timer_split);
                auto stack_size = stack.size();
                if (temp_stack != nullptr) {
                    TIMER_TIC(m_timer_stack_copy);
                    std::copy(std::cbegin(*temp_stack),
                              std::cend(*temp_stack),
                              std::front_inserter(stack));

                    TIMER_TOC(m_timer_stack_copy);
                    forward_fast(ss.str(), index + 1, block, result, stack);
                    // assert(stack.size() == stack_size);
                }
            }
        }
            break;
        case 'T': {
            std::size_t size = shape.getNSamples();

            auto block_offset = calc_offset(offset, block_stride);
            std::unique_ptr<float[]> cache_block;
            std::unique_ptr<float[]> cache_transformed;
            std::unique_ptr<float[]> cache_result;
            auto block_stride_bkp = block_stride;
            auto flat_size_bkp = flat_size;

            const float *block_in;
            float *block_transf;
            float *block_out;

            bool use_cache = (block_offset > 0) && codec_parameters.experimental;

            if (use_cache) {
                TIMER_TIC(m_timer_unique_ptr_alloc);
                cache_block = std::unique_ptr<float[]>(new float[size]);
                cache_transformed = std::unique_ptr<float[]>(new float[size]);
                cache_result = std::unique_ptr<float[]>(new float[size]);
                TIMER_TOC(m_timer_unique_ptr_alloc);

                float *block_in_writable = cache_block.get();
                block_transf = cache_transformed.get();

                TIMER_TIC(m_timer_cache_copy_overhead);
                FOREACH_4D_IDX(i, shape, block_stride)*block_in_writable++ = block[i + block_offset];
                TIMER_TOC(m_timer_cache_copy_overhead);

                block_in = cache_block.get();
                block_out = cache_result.get();
                block_stride = make_stride(shape);
                flat_size = size;
            } else {
                block_in = block;
                block_transf = m_temp_tf_block.get();
                block_out = m_temp_r_block.get();
                if (size > flat_size) {
                    flat_size_bkp = flat_size;
                    flat_size = size;
                }
            }

            for (const auto type : T_CHOICES) {
                std::stringstream ss;

                auto transform_type = static_cast<TransformType>(type);
                ss << prefix << type << suffix;

                TIMER_TIC(m_timer_md_forward_fast);
                md_forward(transform_type, block_in, block_transf, shape);
                TIMER_TOC(m_timer_md_forward_fast);

                auto *ptr_block_transf = block_transf;
                for_each_4d(q_volume, shape, make_stride(shape), [&](auto q){
                    *ptr_block_transf = std::truncf(*ptr_block_transf / q) * q;
                });

                TIMER_TIC(m_timer_md_inverse_fast);
                md_inverse(transform_type, block_transf, block_out, shape);
                TIMER_TOC(m_timer_md_inverse_fast);

                if (use_cache) {
                    float *result = m_temp_r_block.get();
                    float *pout = block_out;
                    TIMER_TIC(m_timer_cache_copy_overhead);
                    FOREACH_4D_IDX(i, shape, block_stride_bkp)
                        result[i + block_offset] = *pout++;
                    TIMER_TOC(m_timer_cache_copy_overhead);
                }

                if (index < descriptor.size() - 1) {
                    forward_fast(ss.str(), index + 1, block, result, stack);
                } else {
                    /* Here we have a full block transformed. */
                    auto old_stride = block_stride;
                    auto old_size = flat_size;
                    block_stride = block_stride_bkp;
                    flat_size = flat_size_bkp;
                    TIMER_TIC(m_timer_rd_cost);
                    calculate_rd_cost(block, ss.str());
                    TIMER_TOC(m_timer_rd_cost);
                    block_stride = old_stride;
                    flat_size = old_size;
                }
            }
            block_stride = block_stride_bkp;
            flat_size = flat_size_bkp;
        }
            break;
    }
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




// EOF