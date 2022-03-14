#ifndef LF_CODEC_ARITHMETICENTROPYENCODER_H
#define LF_CODEC_ARITHMETICENTROPYENCODER_H

#include "EncSymbol.h"
#include "EncoderParameters.h"
#include "EntropyReport.h"
#include "SubpartitionModel.h"

class ArithmeticEntropyEncoder : public EncSymbol {
  public:
    ArithmeticEntropyEncoder(EncoderParameters *parameters, uint bufferSize);
    ~ArithmeticEntropyEncoder();

    void encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, std::string channel);
    void finish_and_write();
    void write_completedBytes();
    uint getTotalBytes() const;

  private:
    // functions
    void open_file(const std::string &filename);

    // variables
    SubpartitionModel *subpartitionModel = nullptr;

    uint totalBytes{0};

    int hypercube = 0;
    std::string ch = "";

    EncoderParameters *parameters;

    std::ofstream outputFile;

    // report
    EntropyReport report;

};

#endif // LF_CODEC_ARITHMETICENTROPYENCODER_H
