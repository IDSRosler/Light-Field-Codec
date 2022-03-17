#ifndef LF_CODEC_ENTROPYREPORT_H
#define LF_CODEC_ENTROPYREPORT_H

#include <fstream>
#include <utility>

struct EntropyReport {

public:
    void openFiles(const std::string& output_path) {
        //this->statisticsFile.open(output_path + "Entropy_Statistics.csv");
        this->subpartitionFile.open(output_path + "reports/tree_csv_report.csv");
        this->treeFlagsFile.open(output_path + "reports/tree_flags.txt");
    }

    void closeFiles() {
        //if (this->statisticsFile.is_open()) this->statisticsFile.close();
        if (this->subpartitionFile.is_open()) this->subpartitionFile.close();
        if (this->treeFlagsFile.is_open()) this->treeFlagsFile.close();
    }

    void setHeaders() {
//        this->statisticsFile << "Hypercube, "
//                                "Channel, "
//                                "Last_Block_Level, "
//                                "Sig_Subpartitions, "
//                                "Non_Sig_Subpartitions, "
//                                "Sig_Coefficients, "
//                                "Non_Sig_Coefficients, "
//                                "One_Coefficients, "
//                                "Two_Coefficients,"
//                                "Gr_Two,"
//                                "Abs_Max_Value,"
//                                "Abs_Mean_Value" << std::endl;

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

/*    void writeStatistics(int last, int sig_sub, int n_sig_sub,
                         int sig_coef, int n_sig_coef, int one, int two, int gr_two,
                         int max, float mean){
        this->statisticsFile <<
                             this->hy  << "," <<
                             this->ch << "," <<
                             last << "," <<
                             sig_sub << "," <<
                             n_sig_sub << "," <<
                             sig_coef << "," <<
                             n_sig_coef << "," <<
                             one << "," <<
                             two << "," <<
                             gr_two << "," <<
                             max << "," <<
                             mean << std::endl;
    }*/

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

private:
    //std::ofstream statisticsFile;
    std::ofstream subpartitionFile;
    std::ofstream treeFlagsFile;
    int hy;
    std::string ch;

};

#endif //LF_CODEC_ENTROPYREPORT_H
