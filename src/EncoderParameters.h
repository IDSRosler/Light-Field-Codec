
#ifndef ENCODERPARAMETERS_H
#define ENCODERPARAMETERS_H


#include <string>

#include "Point4D.h"


class EncoderParameters {
    float qp{1.0};
    std::string path_input{"input/"},
            path_output{"output/"};
    //EDUARDO BEGIN
    std::string prediction{};
    //EDUARDO END

public:
    Point4D dim_LF{0, 0, 0, 0},
            dim_block{0, 0, 0, 0},
            quant_weight_100{1 * 100, 1 * 100, 1 * 100, 1 * 100};
    bool lytro{false};

    bool isLytro() const;

    EncoderParameters() = default;

    void parse_cli(int argc, char **argv);

    float getQp() const;

    const std::string &getPathInput() const;

    const std::string &getPathOutput() const;

    //EDUARDO BEGIN
    const std::string &getPrediction() const;
    //EDUARDO END

    void report();

};


#endif //ENCODERPARAMETERS_H
