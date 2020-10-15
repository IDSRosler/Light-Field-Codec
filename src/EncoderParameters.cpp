
#include "EncoderParameters.h"
#include "Transform.h"

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
        } else if (flag == "-use-transforms") {
            while(argv[++it][0] != '-')
                use_transforms.emplace_back(argv[it]);
            it--;
        } else if (flag == "-transform-min-angular-size") {
            transform_min_angular_size = strtol(argv[++it], nullptr, 10);
        } else if (flag == "-transform-min-spatial-size") {
            transform_min_spatial_size = strtol(argv[++it], nullptr, 10);
        } else if (flag == "-partition-tree-max-depth") {
            Transform::PARTITION_TREE_DEPTH = strtol(argv[++it], nullptr, 10);
        } else if (flag == "-lossless") {
            lossless = true;
        } else if (flag == "-uniform-quantization") {
            uniform_quantization = true;
      //  } else if (flag == "") {
        } else {
            std::cout << "Unused Option: " << argv[it];
            std::cout << "\t" << argv[++it] << std::endl;
        }

    }

    int count = 0;
    std::stringstream ss;
    for (const auto& transform : use_transforms) {
        if (count++ > 0)
            ss << ", ";
        ss << transform;
        Transform::T_CHOICES.push_back(Transform::get_type(transform));
    }
    for (std::size_t n = 0; n <= quadtree_max_inner_nodes; n++)
        Transform::QUADTREE_NODES_COUNT.push_back(4 * n + 1);

    transforms_in_use = ss.str();
    this->dim_LF.updateNSamples();
    this->dim_block.updateNSamples();
    set_global_parameters();
}

float EncoderParameters::getQp() const { return this->qp; }

const std::string &EncoderParameters::getPathInput() const { return this->path_input; }

const std::string &EncoderParameters::getPathOutput() const { return this->path_output; }



void EncoderParameters::report() {
    std::cout <<"[Report]\n";
    display_report(std::cout, "Input path", this->path_input);
    display_report(std::cout, "Output path", this->path_output);
    display_report(std::cout, "LF Dim", this->dim_LF);
    display_report(std::cout, "Block Dim", this->dim_block);
    display_report(std::cout, "Transform", transforms_in_use);
    display_report(std::cout, "Transform Min angular size", transform_min_angular_size);
    display_report(std::cout, "Transform Min spatial size", transform_min_spatial_size);
    display_report(std::cout, "Quadtree max inner nodes", this->quadtree_max_inner_nodes);
    display_report(std::cout, "QP", this->qp);
    display_report(std::cout, "Quantization Weight (*100)", this->quant_weight_100);
    display_report(std::cout, "Lytro", (this->lytro ? "YES" : "NO"));
    display_report(std::cout, "Lambda", this->lambda);
    display_report(std::cout, "Experimental Features", (this->experimental ? "YES" : "NO"));
    std::cout.flush();
}

bool EncoderParameters::isLytro() const { return this->lytro; }
