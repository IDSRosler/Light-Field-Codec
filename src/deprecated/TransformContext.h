#if !defined(__TRANSFORMCONTEXT_H__)
#define __TRANSFORMCONTEXT_H__

#include "Point4D.h"
#include <cstdlib>

/**
 * @brief Abstract class to define the interface to be
 * implemented for all transforms.
 * 
 * @tparam T Numeric type for inner elements.
 */
template <typename T> class TransformContext {
public:
  /** Alias for inner type of elements. */
  using value_type = T; 

  /** Array of sizes for each dimension. */
  const size_t *size; 

  /** Array of strides for each dimension. */
  const size_t *stride;

  /**
   * @brief Performs forward transform on `input`.
   * @param input pointer to memory region to be transformed.
   * @param output pointer to allocated memory to store the result.
   */
  virtual void forward(const T *input, T *output) = 0;
  /**
   * @brief Performs forward transform on `input`.
   * @param input pointer to memory region to be transformed.
   * @param output pointer to allocated memory to store the result.
   * @param size array pointing to dimension sizes of both `input` and `output`
   */
  virtual void forward(const T *input, T *output, const size_t *size) = 0;
  /**
   * @brief Performs inverse transform on `input`.
   * @param input pointer to memory region to be transformed.
   * @param output pointer to allocated memory to store the result.
   */
  virtual void inverse(const T *input, T *output) = 0;
  /**
   * @brief Performs inverse transform on `input`.
   * @param input pointer to memory region to be transformed.
   * @param output pointer to allocated memory to store the result.
   * @param size array pointing to dimension sizes of both `input` and `output`
   */
  virtual void inverse(const T *input, T *output, const size_t *size) = 0;

  /* Virtual destructor */
  virtual ~TransformContext() {}
};

#endif // __TRANSFORMCONTEXT_H__
