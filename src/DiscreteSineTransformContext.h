#if !defined(__DISCRETESINETRANSFORMCONTEXT_H__)
#define __DISCRETESINETRANSFORMCONTEXT_H__

#include <map>
#include <stdexcept>
#include <cmath>
#include "TransformContext.h"
#include "OrthogonalTransformMixin.h"

template <typename T>
class _DiscreteSineTransformContext : public TransformContext<T>
{
protected:
    static std::map<int, const T *> coeff_cache;

    /* Implementation of the helper method to calculate the DCT coefficients */
    static const T *generate_coefficients(const size_t size)
    {
        auto N = size;
        auto *output = new T[N * N];
        auto *p_output = output;
        double s = sqrt(2.0L / (N + 1.0L));
        for (int k = 0; k < N; k++)
        {
            for (int n = 0; n < N; n++)
            {
                double theta = (M_PI * (k + 1.0L) * (n + 1.0L)) / (N + 1.0L);
                *p_output++ = s * sin(theta);
            }
        }
        return output;
    }
};

/** Declaration of static member */
template <typename T>
std::map<int, const T *> _DiscreteSineTransformContext<T>::coeff_cache;

/** Exports the name DiscreteCosineTransformContext with the mixin 
 * OrthogonalTransformContext */
template <typename T>
using DiscreteSineTransformContext = OrthogonalTransformMixin<_DiscreteSineTransformContext<T>>;

#endif // __DISCRETESINETRANSFORMCONTEXT_H__
