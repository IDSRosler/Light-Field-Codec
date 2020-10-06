
#include "EncoderParameters.h"

void EncoderParameters::parse_cli(int argc, char **argv) {
    for (int it = 1; it < argc; ++it) {
        std::string flag = argv[it];

        if (flag == "-blx") {
            this->dim_block.x = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-bly") {
            this->dim_block.y = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-blu") {
            this->dim_block.u = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-blv") {
            this->dim_block.v = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-qx") {
            this->quant_weight_100.x = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qy") {
            this->quant_weight_100.y = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qu") {
            this->quant_weight_100.u = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qv") {
            this->quant_weight_100.v = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-lfx") {
            this->dim_LF.x = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfy") {
            this->dim_LF.y = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfu") {
            this->dim_LF.u = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfv") {
            this->dim_LF.v = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lytro") {
            this->lytro = (short) strtol(argv[++it], nullptr, 10);
        } else if (flag == "-input") {
            this->path_input = argv[++it];
        } else if (flag == "-output") {
            this->path_output = argv[++it];
        } else if (flag == "-qp") {
            this->qp = strtof(argv[++it], nullptr);
        } else {
            std::cout << "Unused Option: " << argv[it];
            std::cout << "\t" << argv[++it] << std::endl;
        }
    }

    this->dim_LF.updateNSamples();
    this->dim_block.updateNSamples();
}

float EncoderParameters::getQp() const {
    return this->qp;
}

const std::string &EncoderParameters::getPathInput() const {
    return this->path_input;
}

const std::string &EncoderParameters::getPathOutput() const {
    return this->path_output;
}

void EncoderParameters::report() {
    std::cout << "################ Encoder Parameters ################" << std::endl;
    std::cout << "Input path: " << this->path_input << std::endl;
    std::cout << "Output path: " << this->path_output << std::endl;
    std::cout << std::endl;
    std::cout << "LF Dim: " << this->dim_LF << std::endl;
    std::cout << "Block Dim: " << this->dim_block << std::endl;
    std::cout << std::endl;
    std::cout << "QP: " << this->qp << std::endl;
    std::cout << "Quantization Weight (*100): " << this->quant_weight_100 << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << std::endl;
}

bool EncoderParameters::isLytro() const {
    return this->lytro;
}

