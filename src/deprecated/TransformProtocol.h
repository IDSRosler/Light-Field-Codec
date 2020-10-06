#if !defined(__TRANSFORM_PROTOCOL_H__)
#define __TRANSFORM_PROTOCOL_H__

#include "BlockCompare.h"
#include "DiscreteCosineTransformContext4D.h"
#include "DiscreteSineTransformContext4D.h"
#include "Point4D.h"
#include "TransformContext.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

template <typename T> class TransformProtocol {
  enum TransformProtocolEnum {
    DCT,
    DST,
    BESTMATCH,
    RANDOM,
  };

  size_t *size;
  size_t *stride;
  TransformProtocolEnum protocol;

  TransformProtocolEnum _dct_fwd(const T *input, T *output, const size_t *size);
  TransformProtocolEnum _dst_fwd(const T *input, T *output, const size_t *size);
  TransformProtocolEnum _bestmatch_fwd(const T *input, T *output,
                                       const size_t *size);

public:
  TransformProtocol(std::string protocol, Point4D &block_size);
  ~TransformProtocol();

  auto forward(const T *input, T *output);
  auto forward(const T *input, T *output, const size_t *size);

  void inverse(const T *input, T *output, const TransformProtocolEnum p);
  void inverse(const T *input, T *output, const TransformProtocolEnum p,
               const size_t *size);
};

template <typename T> TransformProtocol<T>::~TransformProtocol() {
  delete[] this->size;
  delete[] this->stride;
}

template <typename T>
TransformProtocol<T>::TransformProtocol(std::string protocol,
                                        Point4D &block_size) {
  if (protocol == "DCT")
    this->protocol = DCT;
  else if (protocol == "DST")
    this->protocol = DST;
  else if (protocol == "BESTMATCH")
    this->protocol = BESTMATCH;
  else if (protocol == "RANDOM")
    this->protocol = RANDOM;
  else
    throw std::out_of_range("Unknown transform protocol.");
  size = block_size.to_array();
  stride = new size_t[4];
  stride[0] = 1;
  stride[1] = size[0];
  stride[2] = stride[1] * size[2];
  stride[3] = stride[2] * size[3];
}

template <typename T>
typename TransformProtocol<T>::TransformProtocolEnum
TransformProtocol<T>::_dct_fwd(const T *input, T *output, const size_t *size) {
  DiscreteCosineTransformContext4D<T> ctx(this->size, this->stride);
  ctx.forward(input, output, size);
  return DCT;
}

template <typename T>
typename TransformProtocol<T>::TransformProtocolEnum
TransformProtocol<T>::_dst_fwd(const T *input, T *output, const size_t *size) {
  DiscreteSineTransformContext4D<T> ctx(this->size, this->stride);
  ctx.forward(input, output, size);
  return DST;
}

template <typename T>
typename TransformProtocol<T>::TransformProtocolEnum
TransformProtocol<T>::_bestmatch_fwd(const T *input, T *output,
                                     const size_t *size) {
  // TODO: check for border cases. It might not calculate the correct values.
  const size_t FLAT_SIZE = this->stride[3] * this->size[3];
  T dct_out[FLAT_SIZE];
  T dst_out[FLAT_SIZE];

  std::fill(dct_out, dct_out + FLAT_SIZE, 0);
  std::fill(dst_out, dst_out + FLAT_SIZE, 0);

  DiscreteCosineTransformContext4D<T> dct(this->size, this->stride);
  DiscreteSineTransformContext4D<T> dst(this->size, this->stride);
  dst.forward(input, dst_out, size);
  dct.forward(input, dct_out, size);
  if (BlockComparator<T>::compare(dst_out, dct_out, FLAT_SIZE) < 0) {
    std::copy(dst_out, dst_out + FLAT_SIZE, output);
    return DST;
  } else {
    std::copy(dct_out, dct_out + FLAT_SIZE, output);
    return DCT;
  }
}

template <typename T>
auto TransformProtocol<T>::forward(const T *input, T *output) {
  return forward(input, output, this->size);
}

template <typename T>
auto TransformProtocol<T>::forward(const T *input, T *output,
                                   const size_t *size) {
  switch (this->protocol) {
  case DCT:
    return _dct_fwd(input, output, size);
  case DST:
    return _dst_fwd(input, output, size);
  case BESTMATCH:
    return _bestmatch_fwd(input, output, size);
  }
}

template <typename T>
void TransformProtocol<T>::inverse(const T *input, T *output,
                                   const TransformProtocolEnum p) {
  inverse(input, output, p, this->size);
}

template <typename T>
void TransformProtocol<T>::inverse(const T *input, T *output,
                                   const TransformProtocolEnum p,
                                   const size_t *size) {
  switch (p) {
  case DCT: {
    DiscreteCosineTransformContext4D<T> dct(this->size, this->stride);
    dct.inverse(input, output, size);
    break;
  }
  case DST: {
    DiscreteSineTransformContext4D<T> dst(this->size, this->stride);
    dst.inverse(input, output, size);
    break;
  }
  }
}

#endif // __TRANSFORM_PROTOCOL_H__
