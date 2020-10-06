#include "EncBitstreamBuffer.h"

EncBitstreamBuffer::EncBitstreamBuffer(int buffer_size) : buffer_size(buffer_size), byte_pos(0), bits_to_go(8) {
    this->buffer = new Byte[buffer_size];
}

void EncBitstreamBuffer::reset() { this->byte_pos = 0, this->bits_to_go = 8, this->byte_buf = 0; }

EncBitstreamBuffer::~EncBitstreamBuffer() {
//    free(this->buffer);
}

uint EncBitstreamBuffer::getSize() {
    return this->byte_pos * 8 + (8 - this->bits_to_go);
}

EncBitstreamBuffer &EncBitstreamBuffer::operator=(const EncBitstreamBuffer &orig) {
    if (&orig == this) return *this;

    this->buffer_size = orig.buffer_size;
    this->byte_pos = orig.byte_pos;
    this->bits_to_go = orig.bits_to_go;
    this->byte_buf = orig.byte_buf;

    std::copy(orig.buffer, orig.buffer + byte_pos, this->buffer);
    return *this;
}
