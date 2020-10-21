#if !defined(__DISCRETECOSINETRANSFORMCONTEXT_H__)
#define __DISCRETECOSINETRANSFORMCONTEXT_H__

#include "OrthogonalTransformMixin.h"
#include "TransformContext.h"
#include <cmath>
#include <map>
#include <stdexcept>

/**
 * @brief Incomplete class that implements generate_coefficients
 * specific for the Discrete Cosine Transform.
 *
 * @tparam T Numeric type for inner elements.
 */
template <typename T>
class _DiscreteCosineTransformContext : public TransformContext<T> {
protected:
  /** Store previously calculated DCT coefficients for further reuse. */
  static std::map<int, const T *> coeff_cache;

  /**
   * @brief Generate the DCT coefficients matrix.
   * 
   * @param size size of the array to be transformed.
   * @return const T* pointer to an array of `size`^2 elements.
   */
  static const T *generate_coefficients(const size_t size) {
    auto N = size;
    auto *output = new T[N * N];
    auto *p_output = output;
    for (int k = 0; k < N; k++) {
      double s = (k == 0) ? 1 / (double)sqrt(N) : (sqrt(((double)2 / N)));

      for (int n = 0; n < N; n++)
        *p_output++ = s * (double)cos((double)M_PI * (2 * n + 1) * k / (2 * N));
    }
    return output;
  }
};

/** Specialization and declaration of DCT cache. */
template <typename T>
std::map<int, const T *> _DiscreteCosineTransformContext<T>::coeff_cache;

/**
 * @brief The *complete* Discrete Cosine Transform Class implementation.
 * @tparam T Numeric type for inner elements.
 */
template <typename T>
using DiscreteCosineTransformContext =
    OrthogonalTransformMixin<_DiscreteCosineTransformContext<T>>;

#endif // __DISCRETECOSINETRANSFORMCONTEXT_H__
