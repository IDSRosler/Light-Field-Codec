
#ifndef ENCODERPARAMETERS_H
#define ENCODERPARAMETERS_H

#include <string>

#include "Point4D.h"
#include "utils.h"

class EncoderParameters {
    //EDUARDO BEGIN
    std::string prediction{};
    //EDUARDO END

  std::string path_input{"input/"}, path_output{"output/"};

  void set_global_parameters() const;

public:
    bool initialized = false;
    static EncoderParameters parameters;
    float qp{1.0};
    Point4D dim_LF{0, 0, 0, 0},
            dim_block{0, 0, 0, 0},
            quant_weight_100{1 * 100, 1 * 100, 1 * 100, 1 * 100};
    bool lytro{true};
    float lambda = 1;
    float lee_c = 10;
    float lee_ai = 1023;
    float lee_a0 = 1023;
    float lee_bi = 0.04;
    float lee_b0 = 0.01;

    bool show_progress_bar = false;
    bool verbose = false;
    bool experimental = false;
    bool fast = true;

    bool calculate_metrics = true;
    bool display_stages = true;
    bool lossless = false;
    bool uniform_quantization = false;

    bool enable_transforms = true;

    std::vector<std::string> use_transforms;
    std::string transforms_in_use;
    std::size_t quadtree_max_inner_nodes = 1;
    std::size_t transform_min_angular_size = 4;
    std::size_t transform_min_spatial_size = 4;

    bool isLytro() const;

  EncoderParameters() = default;

  void parse_cli(int argc, char *argv[]);

  float getQp() const;

  const std::string &getPathInput() const;

  const std::string &getPathOutput() const;

    //EDUARDO BEGIN
    const std::string &getPrediction() const;
    //EDUARDO END

    void report();

};



#endif // ENCODERPARAMETERS_H
