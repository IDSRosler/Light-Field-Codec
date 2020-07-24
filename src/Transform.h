
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Point4D.h"
#include "Typedef.h"
#include "EncoderParameters.h"
#include "EncBitstreamWriter.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <deque>
#include <memory>

struct Node {
    std::vector<Node *> children;
    int transform_type;

    Node() {
        children.resize(16);
        transform_type = -1;
    };

    ~Node() {
        for (auto &it: children)
            if (it) delete it;
    }

    void set_child(int child_pos, Node *node) { children[child_pos] = node; }

    void set_transform_type(int type) { transform_type = type; }

    void from_string(const char *tree) {
        if (*tree == '\0')
            return;
        else if (tree[1] == '\0') {
            transform_type = (int)(*tree - '0');
        } else {
            std::deque<Node *> queue;
            queue.push_back(this);
            const char *c = tree;
            while (!queue.empty()) {
                auto *node = queue.front();
                queue.pop_front();
                for (int i = 0; i < 16; i++, c++) {
                    auto *child = new Node();
                    if (*c == '0')
                        queue.push_back(child);
                    else
                        child->transform_type = (int)(*c - '0');
                    node->set_child(i, child);
                }
            }
        }
    }

    void to_string(char *buffer) {
        std::deque<Node *> queue;
        queue.push_back(this);
        char *c = buffer;
        while (!queue.empty()) {
            auto *node = queue.front();
            queue.pop_front();
            if (node->transform_type != -1)
                *c++ = (char)(node->transform_type + '0');
            else
                for (int i = 0; i < 16; i++, c++) {
                    auto child = node->children[i];
                    if (child && child->transform_type == -1) {
                        *c = '0';
                        queue.push_back(child);
                    } else if (child) {
                        *c = (char)child->transform_type + '0';
                    }
                }
        }
        *c = '\0';
    }
};

class Transform {
  public:
    enum TransformType {
        NO_TRANSFORM = 0, /* Empty value */
        DCT_II,           /* Discrete Cosine Transform Type II */
        DST_I,            /* Discrete Sine Transform Type II */
        DST_VII,          /* Discrete Sine Transform Type VII */
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

    EncoderParameters codec_parameters;
    bool disable_segmentation = false;

    static void flush_cache();
    static TransformType get_type(std::string transform);

    Transform(Point4D &shape);
    Transform(EncoderParameters &params);

    ~Transform();

    void use_statistics(const std::string filename);
    void set_position(int channel, const Point4D &current_pos);

    std::string forward(TransformType transform, float *input, float *output, Point4D &shape);
    void inverse(TransformType transform, float *input, float *output, Point4D &shape);
    void inverse(const std::string tree, float *input, float *output, Point4D &shape);
    void md_forward(TransformType transform, float *input, float *output, Point4D &shape);
    void md_inverse(TransformType transform, float *input, float *output, Point4D &shape);

  private:
    Point4D shape;
    Point4D stride;
    Point4D position;
    int channel;
    size_t flat_size;
    float *partial_values;
    EncBitstreamWriter *fake_encoder;

    TransformType enforce_transform;

    const std::vector<TransformType> ALL_TRANSFORMS = {DCT_II, DST_I, DST_VII};
    const std::vector<TransformType> HYBRID_LST = {DCT_II, DST_VII, DST_VII, DST_VII};

    static std::map<size_t, float *> cache_dct_ii;
    static std::map<size_t, float *> cache_dst_i;
    static std::map<size_t, float *> cache_dst_vii;

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
    void md_fw_axis(int ax, TransformType type, float *input, float *output, Point4D &shape);
    void md_in_axis(int ax, TransformType type, float *input, float *output, Point4D &shape);

    auto get_block_bpp(float *block, Point4D &shape);
    auto calculate_distortion(float *block, float *result, Point4D &shape);
    auto calculate_tree(float *block, float *result, Point4D &shape, int level);
    void reconstruct_from_tree(Node *root, float *input, float *output, Point4D &shape);
    auto get_quantization_procotol(TransformType transform);
    auto get_transform_vector(TransformType transform);
};

#endif // TRANSFORM_H
