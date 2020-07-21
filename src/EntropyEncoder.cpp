#include "EntropyEncoder.h"

EntropyEncoder::EntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncSymbol(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");
}

EntropyEncoder::~EntropyEncoder() {
    if (this->outputFile.is_open()) this->outputFile.close();
}

void EntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block) {
    this->root = this->tree.CreateRoot(bitstream, dim_block);
    this->tree.CreateTree(this->root, dim_block, {0,0,0,0});
    this->tree.DeleteTree(&this->root);
}

void EntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->outputFile.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;

    this->byte_pos = 0;
}

uint EntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

void EntropyEncoder::open_file(const string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}