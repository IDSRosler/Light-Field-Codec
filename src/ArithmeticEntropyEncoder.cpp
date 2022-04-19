//
// Created by idsrosler on 13/03/2022.
//

#include "ArithmeticEntropyEncoder.h"

ArithmeticEntropyEncoder::ArithmeticEntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncodeSymbol(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");

    // reports files
    this->report.openFiles(this->parameters->getPathOutput());
    this->report.setHeaders();
  }

void ArithmeticEntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, std::string channel) {
    this->hypercube = hypercube_pos;
    this->ch = std::move(channel);

    this->report.setAtt(this->hypercube, this->ch);
    this->report.writeTreeHeader(hypercube_pos, this->ch);

    //make tree partition
    this->treeFlags.clear();
    this->subpartitionModel = new SubpartitionModel(bitstream, dim_block, this->treeFlags, &this->report, hypercube_pos);

    //apply lre in tree header flags
    std::vector<LRE_Struct> lreValues = this->lre->encodeLREVector(this->treeFlags);
    this->encodeLREVector(lreValues);

    this->report.endTreeFlagLine();

    this->subpartitionModel->DeleteTree(); // delete tree
}

void ArithmeticEntropyEncoder::finish_and_write() {
    this->encodeFinalBits();
    this->write_completedBytes();
}

void ArithmeticEntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->outputFile.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;

    this->byte_pos = 0;
}

uint ArithmeticEntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

ArithmeticEntropyEncoder::~ArithmeticEntropyEncoder() {
    this->report.closeFiles();
    this->treeFlags.clear();

    delete this->subpartitionModel;
    delete this->lre;
}

void ArithmeticEntropyEncoder::open_file(const std::string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}
