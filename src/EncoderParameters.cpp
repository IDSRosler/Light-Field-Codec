
#include "EncoderParameters.h"
EncoderParameters EncoderParameters::parameters;

void EncoderParameters::set_global_parameters() const {
  EncoderParameters::parameters = *this;
}

void EncoderParameters::parse_cli(int argc, char *argv[]) {
    this->initialized = true;
    for (int it = 1; it < argc; ++it) {
        std::string flag = argv[it];

        if (flag == "-blx") {
            this->dim_block.x = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-bly") {
            this->dim_block.y = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-blu") {
            this->dim_block.u = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-blv") {
            this->dim_block.v = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-qx") {
            this->quant_weight_100.x = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qy") {
            this->quant_weight_100.y = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qu") {
            this->quant_weight_100.u = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-qv") {
            this->quant_weight_100.v = 100 * strtof(argv[++it], nullptr);
        } else if (flag == "-lfx") {
            this->dim_LF.x = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfy") {
            this->dim_LF.y = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfu") {
            this->dim_LF.u = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lfv") {
            this->dim_LF.v = (short)strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lytro") {
            this->lytro = true;
        } else if (flag == "-input") {
            this->path_input = argv[++it];
        } else if (flag == "-output") {
            this->path_output = argv[++it];
        } else if (flag == "-qp") {
            this->qp = strtof(argv[++it], nullptr);
        } else if (flag == "-lee-c") {
            this->lee_c = strtof(argv[++it], nullptr);
        } else if (flag == "-lee-a0") {
            this->lee_a0 = strtof(argv[++it], nullptr);
        } else if (flag == "-lee-ai") {
            this->lee_ai = strtof(argv[++it], nullptr);
        } else if (flag == "-lee-b0") {
            this->lee_b0 = strtof(argv[++it], nullptr);
        } else if (flag == "-lee-bi") {
            this->lee_bi = strtof(argv[++it], nullptr);
        } else if (flag == "-lambda") {
            this->lambda = strtof(argv[++it], nullptr);
        } else if (flag == "-show-progress-bar") {
            this->show_progress_bar = true;
        } else if (flag == "-verbose") {   
            this->verbose = true;
        } else if (flag == "-experimental"){
            this->experimental = true;
        } else if (flag == "-transform"){
            this->transform = argv[++it];
        } else if (flag == "-enforce-dct-for-non-luma") {
            this->dct_for_non_luma = true;
        } else if (flag == "-max-splits") {
            this->max_splits = (int)strtol(argv[++it], nullptr, 10);
        } else {
            std::cout << "Unused Option: " << argv[it];
            std::cout << "\t" << argv[++it] << std::endl;
        }
    }

    this->dim_LF.updateNSamples();
    this->dim_block.updateNSamples();
    set_global_parameters();
}

float EncoderParameters::getQp() const { return this->qp; }

const std::string &EncoderParameters::getPathInput() const { return this->path_input; }

const std::string &EncoderParameters::getPathOutput() const { return this->path_output; }

void EncoderParameters::report() {
    std::cout <<"[Report]\n";
    display_report(std::cout, "Input path", this->path_input) ;
    display_report(std::cout, "Output path", this->path_output) ;
    display_report(std::cout, "LF Dim", this->dim_LF) ;
    display_report(std::cout, "Block Dim", this->dim_block) ;
    display_report(std::cout, "Transform", this->transform) ;
    display_report(std::cout, "QP", this->qp) ;
    display_report(std::cout, "Quantization Weight (*100)", this->quant_weight_100) ;
    display_report(std::cout, "Lytro", (this->lytro ? "YES" : "NO")) ;
    display_report(std::cout, "Lambda", this->lambda) ;
    display_report(std::cout, "Max Split", this->max_splits) ;
    display_report(std::cout, "Experimental Features", (this->experimental ? "YES" : "NO"))  ;
    std::cout.flush();
}

bool EncoderParameters::isLytro() const { return this->lytro; }
