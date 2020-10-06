#ifndef LF_CODEC_ENCBITSTREAMWRITER_H
#define LF_CODEC_ENCBITSTREAMWRITER_H

#include "EncBitstreamVLC.h"

class EncBitstreamWriter : public EncBitstreamVLC {

public:
    EncBitstreamWriter(EncoderParameters *p, uint bufferSize);

    void writeHeader();

    ~EncBitstreamWriter();

    void write_completedBytes();

    void finish_and_write();

    EncBitstreamWriter &operator+=(const EncBitstreamWriter &temp);

    EncBitstreamWriter &operator+=(const EncBitstreamVLC &temp) override;

    uint getTotalBytes() const;

private:
    std::ofstream fp;

    void open_file(const string &filename);

    uint totalBytes{0};

    EncoderParameters *p;

};

#endif //LF_CODEC_ENCBITSTREAMWRITER_H
