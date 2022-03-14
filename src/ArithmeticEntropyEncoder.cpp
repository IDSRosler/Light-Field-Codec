//
// Created by idsrosler on 13/03/2022.
//

#include "ArithmeticEntropyEncoder.h"
ArithmeticEntropyEncoder::ArithmeticEntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncSymbol(bufferSize){
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

    this->subpartitionModel = new SubpartitionModel(bitstream, dim_block, &this->report);

    this->subpartitionModel->DeleteTree(); // delete tree
}

void ArithmeticEntropyEncoder::finish_and_write() {

}

void ArithmeticEntropyEncoder::write_completedBytes() {

}

uint ArithmeticEntropyEncoder::getTotalBytes() const {

}

ArithmeticEntropyEncoder::~ArithmeticEntropyEncoder() {
    this->report.closeFiles();
}

void ArithmeticEntropyEncoder::open_file(const std::string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}
