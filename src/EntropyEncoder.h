#ifndef LIGHT_FIELD_CODEC_ENTROPYENCODER_H
#define LIGHT_FIELD_CODEC_ENTROPYENCODER_H

#include "EncSymbol.h"
#include "Tree.h"

using namespace std;

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

    Tree tree;
    Node *root = nullptr;
    uint totalBytes{0};
    std::ofstream outputFile;
    EncoderParameters *parameters;
};

#endif //LIGHT_FIELD_CODEC_ENTROPYENCODER_H
