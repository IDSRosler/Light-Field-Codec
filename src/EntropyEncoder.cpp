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
    vector<SyntacticElements> lfbpu_elements; // light field base processing unit

    this->root = this->tree.CreateRoot(bitstream, dim_block);    // create tree
    this->tree.CreateTree(this->root, dim_block, {0,0,0,0});

    this->tree.ComputeLast(last);   // compute last (block level)
    this->encodeLast(last);

    this->tree.ComputeRun(run, last);  // compute run (block level)
    this->encodeRun(run);
    run.clear();

    this->tree.ComputeSyntacticElements(lfbpu_elements, last);  // compute syntactic elements (coefficients level)

    // Todo: encoding of elements here

    lfbpu_elements.clear();
    this->tree.DeleteTree(&this->root); // delete tree
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