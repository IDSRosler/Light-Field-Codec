
#ifndef ENCODERPARAMETERS_H
#define ENCODERPARAMETERS_H

#include <string>

#include "Point4D.h"
#include "utils.h"

class EncoderParameters {

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
    std::string transform;
    float lambda = 1;
    float lee_c = 10;
    float lee_ai = 1023;
    float lee_a0 = 1023;
    float lee_bi = 0.04;
    float lee_b0 = 0.01;
    bool show_progress_bar = false;
    bool verbose = false;
    bool experimental = true;
    bool dct_for_non_luma = false;
    int max_splits = 0;
    bool calculate_metrics = true;
    bool display_stages = true;
    bool calculate_transform_usage = true;
    bool sort_coefficients = true;
    enum {
      ZIG_ZAG_SCAN_ORDER,
      CUSTOM_SCAN_ORDER,
      Z_ORDER_CURVE,
    };
    int scan_order = CUSTOM_SCAN_ORDER;
    bool isLytro() const;

  EncoderParameters() = default;

  void parse_cli(int argc, char *argv[]);

  float getQp() const;

  const std::string &getPathInput() const;

  const std::string &getPathOutput() const;

  void report();
};



#endif // ENCODERPARAMETERS_H
