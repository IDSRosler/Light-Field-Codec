
#ifndef ENCODERPARAMETERS_H
#define ENCODERPARAMETERS_H

#include <string>

#include "Point4D.h"

class EncoderParameters {
  float qp{1.0};
  std::string path_input{"input/"}, path_output{"output/"};

public:
    Point4D dim_LF{0, 0, 0, 0},
            dim_block{0, 0, 0, 0},
            quant_weight_100{1 * 100, 1 * 100, 1 * 100, 1 * 100};
    bool lytro{false};
    std::string transform;
    int flipaxis {0};
    int segments {0};
    float lee_c = 10;
    float lee_ai = 1023;
    float lee_a0 = 1023;
    float lee_bi = 0.04;
    float lee_b0 = 0.01;
    
    bool isLytro() const;

  EncoderParameters() = default;

  void parse_cli(int argc, char *argv[]);

  float getQp() const;

  const std::string &getPathInput() const;

  const std::string &getPathOutput() const;

  void report();
};

#endif // ENCODERPARAMETERS_H
