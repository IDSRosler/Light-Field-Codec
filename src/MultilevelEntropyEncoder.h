#ifndef LF_CODEC_MULTILEVELENTROPYENCODER_H
#define LF_CODEC_MULTILEVELENTROPYENCODER_H

#include "EncodeSymbolsModel.h"
#include "EncoderParameters.h"
#include "EntropyReport.h"
#include "SubpartitionModel.h"
#include "EntropyStructsModel.h"
#include "EntropyLRE.h"

class MultilevelEntropyEncoder: public EncodeSymbolsModel {
  public:
    MultilevelEntropyEncoder(EncoderParameters *parameters, uint bufferSize);
    ~MultilevelEntropyEncoder();

    void encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, std::string channel);
    void finish_and_write();
    void write_completedBytes();
    uint getTotalBytes() const;

  private:
    // functions
    void open_file(const std::string &filename);
    void encodeSyntacticElements(std::queue<Syntactic_Elements> elem);
    void encodeHeaders();
    void initModels();

    // variables
    SubpartitionModel *subpartitionModel = nullptr;

    EntropyLRE *lre = new EntropyLRE();

    std::vector<bool> treeFlags;
    std::queue<Syntactic_Elements> elements;

    uint totalBytes{0};

    int hypercube = 0;
    std::string ch = "";

    EncoderParameters *parameters;

    std::ofstream outputFile;

    // models
    int lastModelU,
        lastModelV,
        sigModel,
        grt1Model,
        grt2Model,
        remModel,
        signModel;

    // report
    EntropyReport report;

};

#endif // LF_CODEC_MULTILEVELENTROPYENCODER_H
