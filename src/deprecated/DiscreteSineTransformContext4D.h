#if !defined(__DISCRETESINETRANSFORMCONTEXT4D_H__)
#define __DISCRETESINETRANSFORMCONTEXT4D_H__

#include "DiscreteSineTransformContext.h"
#include "MultiDimensionalOrthogonalTransformMixin.h"
#include "TransformContext.h"

/**
 * @brief Definition of the 4D Discrete Sine Transform. 
 * @tparam T Numeric type for inner elements.
 */
template <typename T>
using DiscreteSineTransformContext4D =
    MultiDimensionalOrthogonalTransformMixin<TransformContext<T>,
                                             DiscreteSineTransformContext<T>>;

#endif // __DISCRETESINETRANSFORMCONTEXT4D_H__
