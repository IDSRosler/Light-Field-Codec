
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Point4D.h"
#include "Typedef.h"
#include "EncoderParameters.h"
#include "EncBitstreamWriter.h"
#include "Statistics.h"
#include <map>
#include <fstream>
#include <deque>
#include <memory>
#include <stack>
#include <vector>
#include <string>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <functional>

class Transform
{
public:

    enum TransformType
    {
        NO_TRANSFORM = 0, /* Empty value */
        DCT_II = 1,       /* Discrete Cosine Transform Type II */
        DST_I = 2,        /* Discrete Sine Transform Type II */
        DST_VII = 3,      /* Discrete Sine Transform Type VII */
        HYBRID,           /* Applies different transforms across the dimensions */
        MULTI,            /* Pick transform with smallest distortion rate */
    };

    enum
    {
        NO_AXIS = 0,
        AXIS_X = 1,
        AXIS_Y = 2,
        AXIS_U = 4,
        AXIS_V = 8
    };

    inline static const std::string P_CHOICES = "sa";
    inline static const std::vector<TransformType> T_CHOICES = {
        DCT_II,
        // DST_I,
        // DST_VII,
    };
    inline static const std::vector<int> QUADTREE_NODES_COUNT = {1, 5, 9};

    EncoderParameters codec_parameters;
    bool disable_segmentation = false;

    static void flush_cache();
    static TransformType get_type(const std::string &transform);

    Transform(Point4D &shape);
    Transform(EncoderParameters &params);

    ~Transform();

    void use_statistics(const std::string filename);
    void set_position(int channel, const Point4D &current_pos);

    std::pair<std::string, double> forward(const float *block, float *result, const Point4D &shape);
    void forward(const std::string &descriptor, const float *block, float *result, const Point4D &shape);
    void inverse(const std::string &descriptor, const float *block, float *result, const Point4D &shape);
    void md_forward(const TransformType transform, const float *input, float *output, const Point4D &shape);
    void md_inverse(const TransformType transform, const float *input, float *output, const Point4D &shape);

    void md_forward(const TransformType type,
                    const float *input,
                    float *output,
                    const Point4D &_offset,
                    const Point4D &shape);
    void md_inverse(const TransformType type,
                    const float *input,
                    float *output,
                    const Point4D &_offset,
                    const Point4D &shape);

private:
    Point4D block_shape;
    Point4D block_stride;
    Point4D position;

    int channel;
    size_t flat_size;
    float *partial_values;
    float *temp_inverse;
    EncBitstreamWriter *fake_encoder;

    TransformType enforce_transform;
    std::vector<std::string> tree_repr_vector;

    inline static std::map<size_t, float *> cache_dct_ii;
    inline static std::map<size_t, float *> cache_dst_i;
    inline static std::map<size_t, float *> cache_dst_vii;

    

    static float *sd_dct_ii(size_t size);
    static float *sd_dst_i(size_t size);
    static float *sd_dst_vii(size_t size);

    static auto get_transform(TransformType type, bool is_inverse = false);
    static float *get_coefficients(TransformType type, const size_t size);
    static void sd_forward(TransformType type,
                           const float *in,
                           float *out,
                           const size_t offset,
                           const size_t size);
    static void sd_inverse(TransformType type,
                           const float *in,
                           float *out,
                           const size_t offset,
                           const size_t size);
    void md_forward_single_axis(const int ax, TransformType type, const float *input, float *output, const Point4D &shape);
    void md_inverse_single_axis(const int ax, TransformType type, const float *input, float *output, const Point4D &shape);

    template <typename CallbackFunc>
    void forward_fast(
        const std::string &tree_repr,
        std::size_t index,
        const float *block,
        float *result,
        std::deque<std::pair<Point4D, Point4D>> stack,
        CallbackFunc& callback);

    // std::unique_ptr<stack_type> split_coordinate(const char type, const Point4D &_offset, const Point4D &_shape);

    auto get_quantization_procotol(TransformType transform);
    auto get_transform_vector(TransformType transform);
};

#endif // TRANSFORM_H
