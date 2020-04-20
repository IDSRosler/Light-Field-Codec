#if !defined(__TRANSFORMCONTEXT_H__)
#define __TRANSFORMCONTEXT_H__

#include <cstdlib>
#include "Point4D.h"

template <typename T>
class TransformContext
{
public:
    /* Definition of type */
    typedef T value_type;
    const size_t *size;
    const size_t *stride;
    /* Abstract method to be implemented during inheritance */
    virtual void forward(const T *input, T *output) = 0;
    virtual void forward(const T *input, T *output, const size_t *size) = 0;
    /* Abstract method to be implemented during inheritance */
    virtual void inverse(const T *input, T *output) = 0;
    virtual void inverse(const T *input, T *output, const size_t *size) = 0;
    virtual ~TransformContext() {}
};

#endif // __TRANSFORMCONTEXT_H__
