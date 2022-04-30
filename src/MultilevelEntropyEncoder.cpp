//
// Created by idsrosler on 13/03/2022.
//

#include "MultilevelEntropyEncoder.h"

MultilevelEntropyEncoder::MultilevelEntropyEncoder(EncoderParameters *parameters, uint bufferSize) :
    EncodeSymbolsModel(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");

    // reports files
    this->report.openFiles(this->parameters->getPathOutput());
    this->report.setHeaders();
  }

void MultilevelEntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, std::string channel) {
    this->hypercube = hypercube_pos;
    this->ch = std::move(channel);

    this->report.setAtt(this->hypercube, this->ch);
    this->report.writeTreeHeader(hypercube_pos, this->ch);

    //make tree partition
    this->treeFlags.clear();
    while (this->elements.size() > 0){ this->elements.pop(); }
    this->subpartitionModel = new SubpartitionModel(bitstream, dim_block, this->treeFlags, this->elements, &this->report, hypercube_pos);

    //apply lre in tree header flags
    std::vector<LRE_Struct> lreValues = this->lre->encodeLREVector(this->treeFlags);
    this->encodeLREVector(lreValues);

    this->report.endTreeFlagLine();
    this->report.writeSyntactElements(this->elements);

    this->subpartitionModel->DeleteTree(); // delete tree
}

void MultilevelEntropyEncoder::finish_and_write() {
    this->arithDoneEncoding();
    this->write_completedBytes();
}

void MultilevelEntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->outputFile.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;

    this->byte_pos = 0;
}

uint MultilevelEntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

MultilevelEntropyEncoder::~MultilevelEntropyEncoder() {
    this->report.closeFiles();
    this->treeFlags.clear();

    delete this->subpartitionModel;
    delete this->lre;
}

void MultilevelEntropyEncoder::open_file(const std::string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}
