
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
#include <numeric>

class Transform {
public:

    enum TransformType {
        NO_TRANSFORM = 0, /* Empty value */
        DCT_II = 1,       /* Discrete Cosine Transform Type II */
        DST_I = 2,        /* Discrete Sine Transform Type II */
        DST_VII = 3,      /* Discrete Sine Transform Type VII */
        HYBRID,           /* Applies different transforms across the dimensions */
        MULTI,            /* Pick transform with smallest distortion rate */
    };

    enum {
        NO_AXIS = 0,
        AXIS_X = 1,
        AXIS_Y = 2,
        AXIS_U = 4,
        AXIS_V = 8
    };

    inline static const std::string P_CHOICES = "sa";
    inline static std::vector<TransformType> T_CHOICES;
    inline static std::vector<int> QUADTREE_NODES_COUNT;

    EncoderParameters codec_parameters;
    bool disable_segmentation = false;

    static void flush_cache();

    static TransformType get_type(const std::string &transform);

    explicit Transform(Point4D &shape);

    explicit Transform(EncoderParameters &params);

    void set_position(int channel, const Point4D &current_pos);

    std::pair<std::string, double> forward(const float *block, float *result, const Point4D &shape);

    void forward(const std::string &descriptor, const float *block, float *result, const Point4D &shape);

    void inverse(const std::string &descriptor, const float *block, float *result, const Point4D &shape);

    inline void md_forward(TransformType type, const float *input, float *output, const Point4D &shape);

    inline void md_inverse(TransformType type, const float *input, float *output, const Point4D &shape);

    void md_forward(TransformType type,
                    const float *input,
                    float *output,
                    const Point4D &_offset,
                    const Point4D &shape);

    void md_inverse(TransformType type,
                    const float *input,
                    float *output,
                    const Point4D &_offset,
                    const Point4D &shape);

private:
    Point4D block_shape;
    Point4D block_stride;
    Point4D position;

    int channel = 0;
    size_t flat_size = 0;

    double m_min_rd_cost = std::numeric_limits<double>::infinity();
    std::string m_min_descriptor;
    std::unique_ptr<float[]> m_partial_block = nullptr;
    std::unique_ptr<float[]> m_temp_r_block = nullptr;
    std::unique_ptr<float[]> m_temp_tf_block = nullptr;
    std::unique_ptr<int[]> m_temp_lre_block = nullptr;

    std::unique_ptr<EncBitstreamWriter> fake_encoder = nullptr;
    std::unique_ptr<LRE> lre = nullptr;

    TransformType enforce_transform;
    std::vector<std::string> tree_repr_vector;

    inline static std::map<size_t, float *> cache_dct_ii;
    inline static std::map<size_t, float *> cache_dst_i;
    inline static std::map<size_t, float *> cache_dst_vii;


    static float *sd_dct_ii(size_t size);

    static float *sd_dst_i(size_t size);

    static float *sd_dst_vii(size_t size);


    static float *get_coefficients(TransformType type, std::size_t size);

    static void sd_forward(TransformType type,
                           const float *in,
                           float *out,
                           std::size_t offset,
                           std::size_t size);

    static void sd_inverse(TransformType type,
                           const float *in,
                           float *out,
                           std::size_t offset,
                           std::size_t size);

    void md_forward_single_axis(int ax, TransformType type, const float *input, float *output, const Point4D &shape);

    void md_inverse_single_axis(int ax, TransformType type, const float *input, float *output, const Point4D &shape);


    void forward_fast(
            const std::string &tree_repr,
            std::size_t index,
            const float *block,
            float *result,
            std::deque<std::pair<Point4D, Point4D>> stack);

    void calculate_rd_cost(const float *block, const std::string& descriptor);

};


inline void Transform::md_forward(TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &shape)
{
    md_forward(type, input, output, {0, 0, 0, 0}, shape);
}

inline void Transform::md_inverse(TransformType type,
                           const float *input,
                           float *output,
                           const Point4D &shape)
{
    md_inverse(type, input, output, {0, 0, 0, 0}, shape);
}

#endif // TRANSFORM_H
