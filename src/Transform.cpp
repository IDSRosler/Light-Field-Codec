
#include "Transform.h"
#include "EncBitstreamWriter.h"
#include "LRE.h"
#include "Quantization.h"
#include "ScanOrder.h"
#include "utils.h"

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
// #include <execution>
Transform::Transform(Point4D &shape)
{
    this->block_shape = shape;
    this->enforce_transform = TransformType::NO_TRANSFORM;
    block_stride = make_stride(shape);
    flat_size = block_stride.v * shape.v;
    m_partial_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_r_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_tf_block = std::unique_ptr<float[]>(new float[flat_size]);
    m_temp_lre_block = std::unique_ptr<int[]>(new int[flat_size]);
    lre = std::make_unique<LRE>(flat_size == 15 * 15 * 15 * 15);

    for (const auto &N : QUADTREE_NODES_COUNT) {
        for (const auto &tree : generate_full_binary_trees(N)) {
            tree_repr_vector.push_back(tree->repr());
        }
    }
}

Transform::Transform(EncoderParameters &params) : Transform(params.dim_block)
{
    codec_parameters = params;
    fake_encoder = std::make_unique<EncBitstreamWriter>(&codec_parameters, 1048576, true);
}

float *Transform::get_coefficients(Transform::TransformType type,
                                   const size_t size)
{
  float *coeff;
  switch (type)
  {
  case DST_I:
    try
    {
      coeff = cache_dst_i.at(size);
    }
    catch (...)
    {
      coeff = sd_dst_i(size);
      cache_dst_i[size] = coeff;
    }
    break;
  case DST_VII:
    try
    {
      coeff = cache_dst_vii.at(size);
    }
    catch (...)
    {
      coeff = sd_dst_vii(size);
      cache_dst_vii[size] = coeff;
    }
    break;
  case DCT_II:
  default:
    try
    {
      coeff = cache_dct_ii.at(size);
    }
    catch (...)
    {
      coeff = sd_dct_ii(size);
      cache_dct_ii[size] = coeff;
    }
  }
  return coeff;
}

void Transform::sd_forward(Transform::TransformType type, const float *in,
                           float *out, const size_t stride, const size_t size)
{
  float *coeff = get_coefficients(type, size);
    auto pout = out;
    auto pcoeff = coeff;
    for (std::size_t k = 0; k < size; k++, pout += stride) {
        auto pin = in;
        *pout = 0;
        for (std::size_t n = 0; n < size; n++, pin += stride, pcoeff++)
            *pout += *pin * *pcoeff;
    }
}

void Transform::sd_inverse(Transform::TransformType type, const float *in,
                           float *out, const size_t stride, const size_t size)
{
  float *coeff = get_coefficients(type, size);
    auto pout = out;
    for (std::size_t k = 0; k < size; k++, pout += stride) {
        auto pin = in;
        auto pcoeff = coeff + k;
        *pout = 0;
        for (std::size_t n = 0; n < size; n++, pin += stride, pcoeff += size)
            *pout += *pin * *pcoeff;
    }
  
}



float *Transform::sd_dst_vii(size_t size)
{
  auto *output = new float[size * size];
  auto *pout = output;
  for (std::size_t k = 0; k < size; k++)
  {
    double s = (double)2 / (sqrt(2 * size + 1));
    for (std::size_t n = 0; n < size; n++)
      *pout++ = s * sin((double)M_PI * (n + 1.0) * (2.0 * k + 1) / (2.0 * size + 1));
  }
  return output;
}
float *Transform::sd_dst_i(size_t size)
{
  auto *output = new float[size * size];
  auto *pout = output;
  double s = sqrt(2.0 / (size + 1.0));
  for (std::size_t k = 0; k < size; k++)
  {
    for (std::size_t n = 0; n < size; n++)
    {
      double theta = (M_PI * (k + 1.0) * (n + 1.0)) / (size + 1.0);
      *pout++ = s * sin(theta);
    }
  }
  return output;
}

float *Transform::sd_dct_ii(size_t size)
{
  auto *output = new float[size * size];
  float *pout = output;
  for (std::size_t k = 0; k < size; k++)
  {
    double s = (k == 0) ? (1.0 / sqrt(size)) : (sqrt((2.0 / size)));
    for (std::size_t n = 0; n < size; n++)
      *pout++ = s * cos((double)M_PI * ((2.0 * n + 1) * k) / (2.0 * size));
  }
  return output;
}



void Transform::md_forward_single_axis(const int ax,
                                       const TransformType type,
                                       const float *input,
                                       float *output,
                                       const Point4D &shape)
{
  using index_type = Point4D::value_type;
  //using std::execution::par_unseq;
  int axis[3] = {0};
  auto *p = axis;
  for (int x = 0; x < 4; x++)
    if (x != ax)
      *p++ = x;
  // std::vector<index_type> i_indexes(shape[axis[0]]);
  // std::iota(std::begin(i_indexes), std::end(i_indexes), 0);
  // index_type ax_shape = shape[ax];
  // index_type ax_stride = block_stride[ax];
  
  for (index_type k = 0; k < shape[axis[2]]; ++k)
    for (index_type j = 0; j < shape[axis[1]]; ++j)
    {
      for (index_type i = 0; i < shape[axis[0]]; ++i)
      {
        index_type index = k * block_stride[axis[2]] + 
                           j * block_stride[axis[1]] + 
                           i * block_stride[axis[0]];
        sd_forward(type, input + index, output + index, block_stride[ax], shape[ax]);
      }
      // std::for_each(std::begin(i_indexes), std::end(i_indexes),
      //   [&](index_type i) {
      //     index_type index = k * block_stride[axis[2]] + 
      //                        j * block_stride[axis[1]] + 
      //                        i * block_stride[axis[0]];
      //     sd_forward(type, input + index, output + index, ax_stride, ax_shape);
      //   });
    }
      
}

void Transform::md_inverse_single_axis(const int ax,
                                       const TransformType type,
                                       const float *input,
                                       float *output,
                                       const Point4D &shape)
{
  
  using index_type = Point4D::value_type;
  int axis[3] = {0};
  auto *p = axis;
  for (int x = 0; x < 4; x++)
    if (x != ax)
      *p++ = x;

  for (index_type k = 0; k < shape[axis[2]]; ++k)
    for (index_type j = 0; j < shape[axis[1]]; ++j)
      for (index_type i = 0; i < shape[axis[0]]; ++i)
      {
          index_type index = k * block_stride[axis[2]] +
                             j * block_stride[axis[1]] +
                             i * block_stride[axis[0]];
        sd_inverse(type, input + index, output + index, block_stride[ax], shape[ax]);
      }
}



void Transform::md_forward(const TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &offset_,
                           const Point4D &shape)
{
//    int ax_order[] = {3, 2, 1, 0};
//    Point4D adjusted_shape = shape;
//    const float *pin;
//    float *pout;
//    auto offset = calc_offset(offset_, block_stride);
//
//#if LFCODEC_TRANSFORM_CACHE_SMALL_BLOCKS
//    float *wpin;
//    bool should_use_cache = offset > 0;
//    auto old_block_stride = block_stride;
//    std::unique_ptr<float[]> cached_input = nullptr;
//    std::unique_ptr<float[]> cached_output = nullptr;
//#endif
//
//
//    float *intermediary_values = m_partial_block.get();
//
//    for (int ax = 0; ax < 4; ax++)
//        if (offset_[ax] + shape[ax] > block_shape[ax])
//            adjusted_shape[ax] = block_shape[ax] - offset_[ax];
//
//    adjusted_shape.updateNSamples();
//    std::size_t block_size = flat_size;
//
//#if LFCODEC_TRANSFORM_CACHE_SMALL_BLOCKS
//    if (should_use_cache) {
//        block_size = adjusted_shape.getNSamples();
//        cached_input = std::unique_ptr<float[]>(new float[block_size]);
//        cached_output = std::unique_ptr<float[]>(new float[block_size]);
//
//        wpin = cached_input.get();
//        pin = wpin;
//        pout = cached_output.get();
//
//        std::size_t j = 0;
//        FOREACH_4D_IDX(i, adjusted_shape, block_stride) {
//          wpin[j++] = input[i + offset];
//        }
//        block_stride = make_stride(adjusted_shape);
//    } else {
//
//        pin = input;
//        pout = output;
//    }
//#else
//    pin = input + offset;
//    pout = output + offset;
//#endif
//
//    for (auto i : {0, 1, 2, 3}) {
//        md_forward_single_axis(ax_order[i], type, pin, pout, adjusted_shape);
//#if LFCODEC_TRANSFORM_CACHE_SMALL_BLOCKS
//        pin = intermediary_values;
//        if (i < 3)
//            std::copy(pout, pout + block_size, intermediary_values);
//#else
//        pin = intermediary_values + offset;
//        if (i < 3)
//            std::copy(output, output + block_size, intermediary_values);
//#endif
//    }
    int ax_order[] = {3, 2, 1, 0};
    float *partial_values = m_partial_block.get();
    Point4D adjusted_shape = shape;

    for (int ax = 0; ax < 4; ax++)
        if (offset_[ax] + shape[ax] > block_shape[ax])
            adjusted_shape[ax] = block_shape[ax] - offset_[ax];


    int offset = calc_offset(offset_, block_stride);
    const float *pin = input + offset;
    float *pout = output + offset;

    for (auto i : {0, 1, 2, 3})
    {
        md_forward_single_axis(ax_order[i], type, pin, pout, adjusted_shape);
        pin = partial_values + offset;
        pout = output + offset;
        if (i < 3)
            // TODO: Optimize for smaller sizes
            std::copy(output, output + flat_size, partial_values);
    }
#if LFCODEC_TRANSFORM_CACHE_SMALL_BLOCKS
    if (should_use_cache) {
        block_stride = old_block_stride;
        std::size_t j = 0;
        FOREACH_4D_IDX(i, adjusted_shape, block_stride) {
            output[i + offset] = pout[j++];
        }
    }
#endif


#if LFCODEC_USE_QUANTIZATION
    // BUG: code below could be wrong!! volume and pout may have different strides
//  if (codec_parameters.initialized && LFCODEC_USE_QUANTIZATION) {
//    Quantization q(block_shape, codec_parameters);
//    float *volume = q.get_volume(Quantization::HAIYAN);
//    FOREACH_4D_IDX(i, shape, block_stride) {
//      pout[i] = std::trunc(pout[i] / volume[i]);
//    }
//  }
#endif
}

void Transform::md_inverse(const TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &offset_,
                           const Point4D &shape)
{
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

#if LFCODEC_USE_QUANTIZATION
  float block[flat_size];
  if (codec_parameters.initialized && LFCODEC_USE_QUANTIZATION) {
    Quantization q(block_shape, codec_parameters);
    float *volume = q.get_volume(Quantization::HAIYAN);
    
    const float *pin = block + offset;
    FOREACH_4D_IDX(i, shape, block_stride) {
      block[i + offset] = input[i + offset] * volume[i];
    }
  }
#endif  

  for (int i = 0; i < 4; i++)
  {
    md_inverse_single_axis(ax_order[i], type, pin, pout, adjusted_shape);
    pin = intermediary_values + offset;
    pout = output + offset;
    if (i < 3)
    // TODO: Optimize for smaller sizes
      std::copy(output, output + flat_size, intermediary_values);
  }
}

void Transform::flush_cache() {}

Transform::TransformType Transform::get_type(const std::string &transform)
{
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

void Transform::set_position(int _channel, const Point4D &current_pos)
{
  this->channel = _channel;
  this->position = current_pos;
}

std::pair<std::string, double> Transform::forward(const float *block, float *result, const Point4D &_shape)
{
  float *tf_block = m_temp_tf_block.get();
  for (const auto &tree_repr : tree_repr_vector) {
      std::deque<std::pair<Point4D, Point4D>> stack;
      stack.push_front(std::make_pair(_shape, Point4D(0, 0, 0, 0)));
      forward_fast(tree_repr, 0, block, tf_block, stack);
  }

  forward(m_min_descriptor, block, result, _shape);
  return std::make_pair(m_min_descriptor, m_min_rd_cost);
}

void Transform::calculate_rd_cost(const float *block, const std::string& descriptor) {
    std::string min_descriptor;

    int *lre_block = m_temp_lre_block.get();
    float *r_block = m_temp_r_block.get();
    std::fill(lre_block, lre_block + flat_size, 0);
    double sse = std::inner_product(block,
                                    block + flat_size,
                                    r_block,
                                    0.0,
                                    [](auto acc, auto val) { return acc + val; },
                                    [](auto blk, auto rec) { return std::pow(blk - rec, 2); });

    std::transform(r_block, r_block + flat_size, lre_block, std::truncf);
    auto mse = sse / flat_size;
    auto lre_result = lre->encodeLRE(lre_block, flat_size);
    auto lre_size = fake_encoder->write4DBlock(lre_block, flat_size, lre_result);
    fake_encoder->reset();
    double rd_cost = mse + codec_parameters.lambda * lre_size;
    // std::cout << "[log] Pos(x=" << position.x << ","
    //                        "y=" << position.y << ","
    //                        "u=" << position.u << ","
    //                        "v=" << position.v << ","
    //                        "ch=" << channel << ") "
    //                     "Descriptor(text=" << descriptor << ", "
    //                                 "mse=" << mse << ", "
    //                                 "lre_size=" << lre_size << ", "
    //                                 "rd_cost=" << rd_cost << ")\n";
    // std::cout.flush();
    if (rd_cost < m_min_rd_cost) {
        m_min_rd_cost = rd_cost;
        m_min_descriptor = descriptor;
    }
}

void Transform::forward_fast(
        const std::string& descriptor,
        std::size_t index,
        const float *block,
        float *result,
        std::deque<std::pair<Point4D, Point4D>> stack)
{

    std::string prefix(descriptor.substr(0, index > 0 ? index : 0));
    std::string suffix(descriptor.substr(index + 1, descriptor.size() - 1));

    std::unique_ptr<float[]> cache_block;
    std::unique_ptr<float[]> cache_result;

    const auto symbol = descriptor[index];
    auto [shape, offset] = stack.front();
    stack.pop_front();
    switch (symbol)
    {
        case 'P':
        {
            for (const auto type : P_CHOICES)
            {
                std::stringstream ss;
                ss << prefix << type << suffix;
                auto temp_stack = split_coordinate(type, offset, shape);
                if (temp_stack != nullptr) {
                    std::copy(std::cbegin(*temp_stack),
                              std::cend(*temp_stack),
                              std::front_inserter(stack));

                    forward_fast(ss.str(), index + 1, block, result, stack);
                }
            }
        }
            break;
        case 'T':
        {

            for (const auto type : T_CHOICES)
            {
                std::stringstream ss;
                auto transform_type = static_cast<TransformType>(type);
                ss << prefix << type << suffix;
                md_forward(transform_type, block, result, offset, shape);
                md_inverse(transform_type, result, m_temp_r_block.get(), offset, shape);

                if (index < descriptor.size() - 1) {
                    forward_fast(ss.str(), index + 1, block, result, stack);
                } else {
                    /* Here we have a full block transformed. */
                    calculate_rd_cost(block, ss.str());
                }
            }
        }
            break;
    }
}



void Transform::forward(const std::string& descriptor, const float *block, float *result, const Point4D &_shape)
{
  for (const auto &desc: parse_descriptor(descriptor, _shape, true))
  {
    const auto &[offset, shape, type] = desc;
    md_forward(static_cast<TransformType>(type), block, result, offset, shape);
  }
}

void Transform::inverse(const std::string& descriptor, const float *block, float *result, const Point4D &_shape)
{
  for (const auto &desc: parse_descriptor(descriptor, _shape, true))
  {
    const auto &[offset, shape, type] = desc;
    md_inverse(static_cast<TransformType>(type), block, result, offset, shape);
  }
}



// EOF