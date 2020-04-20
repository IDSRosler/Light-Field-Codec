#if !defined(__ORTHOGONALTRANSFORMMixin_H__)
#define __ORTHOGONALTRANSFORMMixin_H__

/** Orthogonal Transforms follow these two steps:
 * 1) Calculate its coefficients matrix
 * 2) Multiplication between input array and coefficients matrix
 * 
 * To provide the same interface to execute the same steps, this
 * mixin should inherit from the any class that defines only the
 * method to calculate the coefficients.
 * The necessary methods implementation of .forward() and .inverse()
 * needed for TransformContext<T> are done here. 
 * 
 * To use this class, you need to:
 * 1) Inherit from TransformContext
 * 2) Declare a static property std::map<int, const T*> coeff_cache
 * 3) Use as following:  OrthogonalTransformContext<YourType<T>> *ctx;
 */
template <typename BASE, typename T = typename BASE::value_type>
class OrthogonalTransformMixin : public BASE
{
private:
    /** Method to populate BASE::coeff_cache with the coefficients to be used
     * throuhout the application. Ii assumes BASE::generate_coefficients from
     * the base class. */
    static const T *get_coefficients(const size_t size)
    {
        try
        {
            // Just returns whatever is stored at coeff_cache[size].
            return BASE::coeff_cache.at(size);
        }
        catch (const std::out_of_range &e)
        {
            // Expected exception for cases where coeff_cache[size] is invalid.
            const T *coeff = BASE::generate_coefficients(size);
            BASE::coeff_cache[size] = coeff;
            return coeff;
        }
    }

    /** Method that multiplies the input array with the coefficients matrix.
     * It allows for the transpose run of the coefficients matrix to calculate
     * the inverse transform as well. */
    void ortho_transform(const T *coefficients,
                         const T *array,
                         T *result,
                         const bool transpose_matrix)
    {
        auto *size = BASE::size;
        auto *stride = BASE::stride;
        size_t index_coeff;
        T sum;
        for (int k = 0; k < *size; k++)
        {
            sum = 0;
            for (int n = 0; n < *size; n++)
            {
                index_coeff = transpose_matrix ? n * *size + k : k * *size + n;
                sum += array[n * *stride] * coefficients[index_coeff];
            }
            result[k * *stride] = sum;
        }
    }

    const size_t _size;
    const size_t _stride;

public:
    OrthogonalTransformMixin(const size_t size_, const size_t stride_ = 1)
        : _size(size_), _stride(stride_)
    {
        // Update parent pointers to correct values.
        BASE::size = &_size;
        BASE::stride = &_stride;
    }

    OrthogonalTransformMixin(const size_t *size_, const size_t *stride_)
    {
        OrthogonalTransformMixin(*size_, *stride_);
    }

    void forward(const T *input, T *output)
    {
        forward(input, output, BASE::size);
    }
    void forward(const T *input, T *output, const size_t *size)
    {
        auto *coeff = get_coefficients(*size);
        ortho_transform(coeff, input, output, /* dont transpose */ false);
    }

    void inverse(const T *input, T *output)
    {
        inverse(input, output, BASE::size);
    }
    void inverse(const T *input, T *output, const size_t *size)
    {
        auto *coeff = get_coefficients(*size);
        ortho_transform(coeff, input, output, /* transpose */ true);
    }
    // Silently drops every entry in coeff_cache.
    static void flush_coeff()
    {
        for (const auto &kv : BASE::coeff_cache)
            delete[] kv.second;
        BASE::coeff_cache.clear();
    }
};

#endif // __ORTHOGONALTRANSFORMMixin_H__
