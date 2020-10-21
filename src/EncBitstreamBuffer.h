#ifndef LF_CODEC_ENCBITSTREAMBUFFER_H
#define LF_CODEC_ENCBITSTREAMBUFFER_H

typedef unsigned char Byte;
typedef unsigned int uint;

#include <algorithm>

class EncBitstreamBuffer { // bit_stream_enc: JM
public:
    int buffer_size; //!< Buffer size
    uint byte_pos; //!< current position in bitstream;
    uint bits_to_go; //!< current bitcounter

    Byte byte_buf{}; //!< current buffer for last written byte
    Byte *buffer; //!< actual buffer for written bytes

    explicit EncBitstreamBuffer(int buffer_size = 50);

    EncBitstreamBuffer &operator=(const EncBitstreamBuffer &orig);

    ~EncBitstreamBuffer();

    virtual void reset();

    virtual uint getSize();
};


#endif //LF_CODEC_ENCBITSTREAMBUFFER_H
