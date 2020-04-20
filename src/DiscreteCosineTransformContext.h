#if !defined(__DISCRETECOSINETRANSFORMCONTEXT_H__)
#define __DISCRETECOSINETRANSFORMCONTEXT_H__

#include <map>
#include <stdexcept>
#include <cmath>
#include "TransformContext.h"
#include "OrthogonalTransformMixin.h"

template <typename T>
class _DiscreteCosineTransformContext : public TransformContext<T>
{
protected:
    static std::map<int, const T *> coeff_cache;

    /* Implementation of the helper method to calculate the DCT coefficients */
    static const T *generate_coefficients(const size_t size)
    {
        auto N = size;
        auto *output = new T[N * N];
        auto *p_output = output;
        for (int k = 0; k < N; k++)
        {
            double s = (k == 0) ? 1 / (double)sqrt(N) : (sqrt(((double)2 / N)));

            for (int n = 0; n < N; n++)
                *p_output++ = s * (double)cos((double)M_PI * (2 * n + 1) * k / (2 * N));
        }
        return output;
    }
};

/** Declaration of static member */
template <typename T>
std::map<int, const T *> _DiscreteCosineTransformContext<T>::coeff_cache;

/** Exports the name DiscreteCosineTransformContext with the mixin 
 * OrthogonalTransformContext */
template <typename T>
using DiscreteCosineTransformContext = OrthogonalTransformMixin<_DiscreteCosineTransformContext<T>>;

#endif // __DISCRETECOSINETRANSFORMCONTEXT_H__
