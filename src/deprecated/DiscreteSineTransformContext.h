#if !defined(__DISCRETESINETRANSFORMCONTEXT_H__)
#define __DISCRETESINETRANSFORMCONTEXT_H__

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
class _DiscreteSineTransformContext : public TransformContext<T> {
protected:
/** Store previously calculated DST coefficients for further reuse. */
  static std::map<int, const T *> coeff_cache;

  /**
   * @brief Generate the DST coefficients matrix.
   * 
   * @param size size of the array to be transformed.
   * @return const T* pointer to an array of `size`^2 elements.
   */
  static const T *generate_coefficients(const size_t size) {
    auto N = size;
    auto *output = new T[N * N];
    auto *p_output = output;
    double s = sqrt(2.0L / (N + 1.0L));
    for (int k = 0; k < N; k++) {
      for (int n = 0; n < N; n++) {
        double theta = (M_PI * (k + 1.0L) * (n + 1.0L)) / (N + 1.0L);
        *p_output++ = s * sin(theta);
      }
    }
    return output;
  }
};

/** Specialization and declaration of DST cache. */
template <typename T>
std::map<int, const T *> _DiscreteSineTransformContext<T>::coeff_cache;

/**
 * @brief The *complete* Discrete Sine Transform Class implementation.
 * @tparam T Numeric type for inner elements.
 */
template <typename T>
using DiscreteSineTransformContext =
    OrthogonalTransformMixin<_DiscreteSineTransformContext<T>>;

#endif // __DISCRETESINETRANSFORMCONTEXT_H__
