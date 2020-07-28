#include "EntropyEncoder.h"

EntropyEncoder::EntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncSymbol(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");
}

EntropyEncoder::~EntropyEncoder() {
}

void EntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block) {
    int last;
    vector<int> run;
    this->root = this->tree.CreateRoot(bitstream, dim_block);
    this->tree.CreateTree(this->root, dim_block, {0,0,0,0});
    last  = this->tree.ComputeLast();
    this->encodeLast(last);
    run = this->tree.ComputeRun(last);
    this->encodeRun(run);
    run.clear();
    this->tree.DeleteTree(&this->root);
}

void EntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->outputFile.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;

    this->byte_pos = 0;
}

void EntropyEncoder::finish_and_write() {
    this->writeSyntaxElement(0, this->bits_to_go);
    this->write_completedBytes();
    if (this->outputFile.is_open()) this->outputFile.close();
}

uint EntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

void EntropyEncoder::open_file(const string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}