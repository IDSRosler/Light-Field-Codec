#if !defined(__MULTIDIMORTHOGONALTRANSFORMMIXIN_H__)
#define __MULTIDIMORTHOGONALTRANSFORMMIXIN_H__

#define FORWARD 1
#define INVERSE 2

template <typename BASE, typename CTX, typename T = typename BASE::value_type>
class MultiDimensionalOrthogonalTransformMixin
    : public BASE
{
private:
    bool should_delete_pointers = false;
    size_t direction = 0;
    TransformContext<T> *ctx;

    /**
     * Method used to apply some transformation across a single 
     * dimension given by a parameter `axis`
     * As long as `size` is consistent, all strides can be
     * calculated correctly. 
     */
    void ortho_transform_multidim_step(
        const T *input,
        T *output,
        const size_t *size,
        size_t axis)
    {
        // dims are all other dimensions but the one given to the function.
        // The transform applied on the X axis, for instance, will loop
        // through Y, U and V.
        // This array holds, then, all the other dimensions.
        uint dims[3];

        // populate dims
        for (int j = 0, i = 0; i < 4; i++)
        {
            if (i != axis)
                dims[j++] = i;
        }

        // Iterating accross all other axis.
        // Outter loop holds bigger stride steps.
        for (uint d3 = 0; d3 < size[dims[2]]; ++d3)
        {
            for (uint d2 = 0; d2 < size[dims[1]]; ++d2)
            {
                for (uint d1 = 0; d1 < size[dims[0]]; ++d1)
                {
                    uint index = d1 * BASE::stride[dims[0]] + d2 * BASE::stride[dims[1]] + d3 * BASE::stride[dims[2]];
                    switch (direction)
                    {
                    case FORWARD:
                        ctx->forward(input + index, output + index);
                        break;
                    case INVERSE:
                        ctx->inverse(input + index, output + index);
                        break;
                    default:
                        break;
                    }
                }
            }
        }
    }

    void ortho_transform_multidim(
        const T *input,
        T *output,
        const size_t *size)
    {
        int start_at;
        int end_at;
        int step;

        auto *stride = BASE::stride;
        size_t FULL_LENGTH = stride[3] * size[3];
        T *partial_result = new T[FULL_LENGTH];

        switch (direction)
        {
        case FORWARD:
            start_at = 0;
            end_at = 4;
            step = 1;
            break;
        case INVERSE:
            start_at = 3;
            end_at = -1;
            step = -1;
        default:
            break;
        }

        // For the first round, the "last round output" is simply
        // the input values to be calculated.
        T *pout = (T *)input;

        for (size_t dim = start_at; dim != end_at; dim += step)
        {
            ctx = new CTX(size[dim], stride[dim]);

            // Copy all values from last round into`partial_result`.
            std::copy(pout, pout + FULL_LENGTH, partial_result);

            ortho_transform_multidim_step(partial_result, output, size, dim);

            delete ctx;
            pout = output;
        }
        delete[] partial_result;
    }

public:
    ~MultiDimensionalOrthogonalTransformMixin()
    {
        if (should_delete_pointers)
        {
            delete[] BASE::size;
            delete[] BASE::stride;
        }
    }
    MultiDimensionalOrthogonalTransformMixin(size_t *size_, size_t *stride_)
    {
        BASE::size = size_;
        BASE::stride = stride_;
    }
    MultiDimensionalOrthogonalTransformMixin(const Point4D &size_, const Point4D &stride_)
    {
        BASE::size = size_.to_array();
        BASE::stride = stride_.to_array();
        should_delete_pointers = true;
    }

    void forward(const T *input, T *output)
    {
        forward(input, output, BASE::size);
    }

    void forward(const T *input, T *output, const size_t *size)
    {
        direction = FORWARD;
        ortho_transform_multidim(input, output, size);
    }

    void inverse(const T *input, T *output)
    {
        inverse(input, output, BASE::size);
    }
    void inverse(const T *input, T *output, const size_t *size)
    {
        direction = INVERSE;
        ortho_transform_multidim(input, output, size);
    }
};

#endif // __MULTIDIMORTHOGONALTRANSFORMMIXIN_H__
