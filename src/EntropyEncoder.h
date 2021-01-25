#ifndef LIGHT_FIELD_CODEC_ENTROPYENCODER_H
#define LIGHT_FIELD_CODEC_ENTROPYENCODER_H

#include "EncSymbol.h"
#include "ArithmeticEncoder.h"
#include "Tree.h"
#include "EntropyReport.h"

using namespace std;

struct ElementsFrequency {
    vector<int> sig;
    vector<int> gr_one;
    vector<int> gr_two;
    vector<int> sign;

    void setFrequency(int sig0, int sig1, int gr_one0, int gr_one1, int gr_two0, int gr_two1, int sign0, int sign1){
        this->sig.push_back(sig0); //95
        this->sig.push_back(sig1); //5

        this->gr_one.push_back(gr_one0); //45
        this->gr_one.push_back(gr_one1); //55

        this->gr_two.push_back(gr_two0); //25
        this->gr_two.push_back(gr_two1); //75

        this->sign.push_back(sign0); //50
        this->sign.push_back(sign1); //50
    }
};

class EntropyEncoder : public EncSymbol{
public:
    EntropyEncoder(EncoderParameters *parameters, uint bufferSize);
    ~EntropyEncoder();

    void encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, string channel);
    void finish_and_write();
    void write_completedBytes();
    uint getTotalBytes() const;

private:
    void open_file(const string &filename);

    void EncodeSyntacticElements(vector<SyntacticElements> &lfbpu);
    void encodeSyntacticElements(vector<SyntacticElements> &lfbpu);

    void ComputeFrequency(vector<SyntacticElements> &lfbpu, ElementsFrequency& freq);
    void ComputeFrequency(vector<int> &freq);

    int encodeSymbol(int code, int model, const std::string& type);

    Tree tree;
    Node *root = nullptr;

    vector<int> e_buffer;

    uint totalBytes{0};

    std::ofstream outputFile;
    /*std::ofstream statistics_file;
    std::ofstream bitrate_steps;*/

    int hypercube,
        last,
        sig_sub,
        n_sig_sub,
        sig_coeff,
        n_sig_coeff,
        one,
        two,
        gr_two,
        max_value;

    float mean_value;

    int model_8bits;
    int model_exp;
    int model_sig;
    int model_grone;
    int model_grtwo;
    int model_sign;

    string ch;

    std::ofstream freqFile;

    EntropyReport report;

    EncoderParameters *parameters;

    ArithmeticEncoder arith_encoder = ArithmeticEncoder(*this->buffer, this->bits_to_go, this->byte_buf, this->byte_pos);
};

#endif //LIGHT_FIELD_CODEC_ENTROPYENCODER_H
