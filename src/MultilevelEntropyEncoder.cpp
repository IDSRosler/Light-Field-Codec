//
// Created by idsrosler on 13/03/2022.
//

#include "MultilevelEntropyEncoder.h"

MultilevelEntropyEncoder::MultilevelEntropyEncoder(EncoderParameters *parameters, uint bufferSize) :
    EncodeSymbolsModel(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");

    //  add models
    this->lastModelU = this->arithAddModel();
    this->arithStartModel(this->parameters->dim_block.u, this->lastModelU);
    this->lastModelV = this->arithAddModel();
    this->arithStartModel(this->parameters->dim_block.v, this->lastModelV);

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

    this->encodeSyntacticElements(this->elements);

    //restart models
    //this->arithRestartModel(this->lastModelU);
    //this->arithRestartModel(this->lastModelV);

    this->report.endTreeFlagLine();
    this->report.writeSyntactElements(this->elements);

    this->subpartitionModel->DeleteTree(); // delete tree
}

void MultilevelEntropyEncoder::encodeSyntacticElements(std::queue<Syntactic_Elements> elem) {
    Syntactic_Elements blockElements;
    while (!elem.empty()){
        blockElements = elem.front();
        // encode last u
        this->arithEncodeSymbol(blockElements.last_sig_coeff_u, this->lastModelU);
        this->arithUpdateModel(blockElements.last_sig_coeff_u, this->lastModelU);
        //encode last v
        this->arithEncodeSymbol(blockElements.last_sig_coeff_v, this->lastModelV);
        this->arithUpdateModel(blockElements.last_sig_coeff_v, this->lastModelV);
        elem.pop();
    }
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
