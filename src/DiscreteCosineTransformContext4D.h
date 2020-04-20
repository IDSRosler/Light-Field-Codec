#if !defined(__DISCRETECOSINETRANSFORMCONTEXT4D_H__)
#define __DISCRETECOSINETRANSFORMCONTEXT4D_H__

#include "TransformContext.h"
#include "MultiDimensionalOrthogonalTransformMixin.h"
#include "DiscreteCosineTransformContext.h"

template <typename T>
using DiscreteCosineTransformContext4D = 
    MultiDimensionalOrthogonalTransformMixin<TransformContext<T>,
                                            DiscreteCosineTransformContext<T>>;

#endif // __DISCRETECOSINETRANSFORMCONTEXT4D_H__
