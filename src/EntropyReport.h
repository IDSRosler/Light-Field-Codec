#ifndef LF_CODEC_ENTROPYREPORT_H
#define LF_CODEC_ENTROPYREPORT_H

#include <fstream>
#include <utility>

#include "ArithmeticStructures.h"

struct EntropyReport {

public:
    void openFiles(const std::string& output_path) {
        this->subpartitionFile.open(output_path + "reports_tree_csv_report.csv");
        this->treeFlagsFile.open(output_path + "reports_tree_flags.txt");
        this->block4dFile.open(output_path + "reports_block_uv_level_4.csv");
        this->syntacticsElementFile.open(output_path + "reports_syntactic_elements.csv");
    }

    void closeFiles() {
        if (this->subpartitionFile.is_open()) this->subpartitionFile.close();
        if (this->treeFlagsFile.is_open()) this->treeFlagsFile.close();
        if (this->block4dFile.is_open()) this->block4dFile.close();
        if (this->syntacticsElementFile.is_open()) this->syntacticsElementFile.close();
    }

    void setHeaders() {

        this->subpartitionFile << "Hypercube, "
                                  "Channel, "
                                  "Level, "
                                  "Start_Position_x, "
                                  "Start_Position_y, "
                                  "End_Position_x, "
                                  "End_Position_y, "
                                  "Hy_size, "
                                  "N_zero, "
                                  "N_one, "
                                  "N_two, "
                                  "N_gr_two, "
                                  "Abs_max_value, "
                                  "Abs_mean_value, "
                                  "Is_sig" << std::endl;

        this->treeFlagsFile << "(Hypercube - Channel) -> tree flags \n" << std::endl;

        this->syntacticsElementFile << "Hypercube, "
                                  "Channel, "
                                  "Block_index, "
                                  "Last_u, "
                                  "Last_v, "
                                  "N_Significant_values, "
                                  "N_Non_Significant_values, "
                                  "N_abs_level_greater1, "
                                  "N_abs_level_greater2, "
                                  "N_positive_values, "
                                  "N_negative_values" << std::endl;

    }

    void writeSyntactElements(std::queue<Syntactic_Elements> elements){
        int sig=0;
        int n_sig=0;
        int grt1=0;
        int grt2=0;
        int pos=0;
        int neg=0;
        int index=0;
        Syntactic_Elements elem;
        while (!elements.empty()) {
            elem = elements.front();

            while (!elem.sig_coeff_flag.empty()){
                if (elem.sig_coeff_flag.front()){
                    ++sig;
                }else{
                    ++n_sig;
                }
                elem.sig_coeff_flag.pop();
            }
            while (!elem.coeff_abs_level_greater1_flag.empty()){
                if (elem.coeff_abs_level_greater1_flag.front()){
                    ++grt1;
                }
                elem.coeff_abs_level_greater1_flag.pop();
            }
            while (!elem.coeff_abs_level_greater2_flag.empty()){
                if (elem.coeff_abs_level_greater2_flag.front()){
                    ++grt2;
                }
                elem.coeff_abs_level_greater2_flag.pop();
            }
            while (!elem.coeff_sign_flag.empty()){
                if (elem.coeff_sign_flag.front()){
                    ++neg;
                }else{
                    ++pos;
                }
                elem.coeff_sign_flag.pop();
            }

            this->syntacticsElementFile <<
                this->hy << "," <<
                this->ch << "," <<
                index << "," <<
                elem.last_sig_coeff_u << "," <<
                elem.last_sig_coeff_v << "," <<
                sig << "," <<
                n_sig << "," <<
                grt1 << "," <<
                grt2 << "," <<
                pos << "," <<
                neg << std::endl;

            sig=0;
            n_sig=0;
            grt1=0;
            grt2=0;
            pos=0;
            neg=0;
            elements.pop();

            ++index;
        }
    }

    void writeSubpartitions(int level, int px, int py, int pu,
                            int pv, int zero, int one, int two, int gr_two,
                            int max, float mean, bool is_sig, int hy_size){
        this->subpartitionFile <<
        this->hy << "," <<
        this->ch << "," <<
        level << "," <<
        px << "," <<
        py << "," <<
        pu << "," <<
        pv << "," <<
        hy_size << "," <<
        zero << "," <<
        one << "," <<
        two << "," <<
        gr_two << "," <<
        max << "," <<
        mean << "," <<
        is_sig << std::endl;
    }

    void writeTreeHeader(int hypercube, std::string ch){
        this->treeFlagsFile << "( " + std::to_string(hypercube) + " - " + ch + " ) -> ";
    }

    void writeTreeFlag(int flag){
        this->treeFlagsFile << std::to_string(flag);
    }

    void endTreeFlagLine(){
        this->treeFlagsFile << std::endl;
    }

    void setAtt(int hyp, std::string channel){
        this->hy = hyp;
        this->ch = std::move(channel);
    }

    void writeblock(int ****data, NodeBlock *node){

        for (int it_v = 0; it_v < node->dimention.v; ++it_v) {
            for (int it_u = 0; it_u < node->dimention.u; ++it_u) {
                this->block4dFile << data[node->start_index.x][node->start_index.y][it_u][it_v] << ",";
            }
            this->block4dFile << std::endl;
        }

        this->block4dFile << std::endl;
    }

private:
    std::ofstream subpartitionFile;
    std::ofstream treeFlagsFile;
    std::ofstream block4dFile;
    std::ofstream syntacticsElementFile;

    int hy;
    std::string ch;

};

#endif //LF_CODEC_ENTROPYREPORT_H
