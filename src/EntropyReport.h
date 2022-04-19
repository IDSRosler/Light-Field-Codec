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
        this->block4d.open(output_path + "reports_block_uv_level_4.csv");
    }

    void closeFiles() {
        if (this->subpartitionFile.is_open()) this->subpartitionFile.close();
        if (this->treeFlagsFile.is_open()) this->treeFlagsFile.close();
        if (this->block4d.is_open()) this->treeFlagsFile.close();
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

        for (int it_y = node->start_index.y; it_y < node->end_index.y; ++it_y) {
            for (int it_x = node->start_index.x; it_x < node->end_index.x; ++it_x) {
                for (int it_v = 0; it_v < node->dimention.v; ++it_v) {
                    for (int it_u = 0; it_u < node->dimention.u; ++it_u) {
                        this->block4d << data[it_x][it_y][it_u][it_v] << ",";
                    }
                    this->block4d << std::endl;
                }
            }
        }
        this->block4d << std::endl;
    }

private:
    std::ofstream subpartitionFile;
    std::ofstream treeFlagsFile;
    std::ofstream block4d;

    int hy;
    std::string ch;

};

#endif //LF_CODEC_ENTROPYREPORT_H
