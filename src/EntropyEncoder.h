#ifndef LIGHT_FIELD_CODEC_ENTROPYENCODER_H
#define LIGHT_FIELD_CODEC_ENTROPYENCODER_H

#include "EncSymbol.h"
#include "ArithmeticEncoder.h"
#include "Tree.h"

using namespace std;

struct ElementsFrequency {
    vector<int> sig;
    vector<int> gr_one;
    vector<int> gr_two;
    vector<int> sign;

    void setFrequency(){
        this->sig.push_back(95); //95
        this->sig.push_back(5); //5

        this->gr_one.push_back(45); //45
        this->gr_one.push_back(55); //55

        this->gr_two.push_back(25); //25
        this->gr_two.push_back(75); //75

        this->sign.push_back(50); //50
        this->sign.push_back(50); //50
    }
};

class EntropyEncoder : public EncSymbol{
public:
    EntropyEncoder(EncoderParameters *parameters, uint bufferSize);
    ~EntropyEncoder();

    void encodeHypercube(int *bitstream, const Point4D &dim_block);
    void finish_and_write();
    void write_completedBytes();
    uint getTotalBytes() const;

private:
    void open_file(const string &filename);
    void EncodeSyntacticElements(vector<SyntacticElements> lfbpu);

    void ComputeFrequency(vector<SyntacticElements> lfbpu);

    Tree tree;
    Node *root = nullptr;

    uint totalBytes{0};

    std::ofstream outputFile;

    std::ofstream freqFile;

    EncoderParameters *parameters;

    ArithmeticEncoder arith_encoder = ArithmeticEncoder(*this->buffer, this->bits_to_go, this->byte_buf, this->byte_pos);
};

#endif //LIGHT_FIELD_CODEC_ENTROPYENCODER_H
