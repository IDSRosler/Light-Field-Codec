#include "EncBitstreamWriter.h"

EncBitstreamWriter::EncBitstreamWriter(EncoderParameters *p, uint bufferSize) : EncBitstreamVLC(bufferSize, p) {
    this->p = p;
    open_file(this->p->getPathOutput() + "binary.bin");
}

void EncBitstreamWriter::writeHeader() {
    this->writeGolomb_ui(this->p->dim_LF.x);
    this->writeGolomb_ui(this->p->dim_LF.y);
    this->writeGolomb_ui(this->p->dim_LF.u);
    this->writeGolomb_ui(this->p->dim_LF.v);
    this->writeGolomb_ui(this->p->dim_block.x);
    this->writeGolomb_ui(this->p->dim_block.y);
    this->writeGolomb_ui(this->p->dim_block.u);
    this->writeGolomb_ui(this->p->dim_block.v);
    this->writeGolomb_ui(this->p->quant_weight_100.u);
    this->writeGolomb_ui(this->p->quant_weight_100.v);
    this->writeGolomb_ui(this->p->quant_weight_100.x);
    this->writeGolomb_ui(this->p->quant_weight_100.y);
    this->writeGolomb_ui(floor(this->p->getQp() * 100));
}

void EncBitstreamWriter::open_file(const std::string &filename) {
    fp.open(filename, std::ios::binary);
    assert(fp.is_open());
}

void EncBitstreamWriter::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->fp.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;

    this->byte_pos = 0;
}

void EncBitstreamWriter::finish_and_write() {

    this->writeSyntaxElement(0, this->bits_to_go);

    this->write_completedBytes();
}

EncBitstreamWriter::~EncBitstreamWriter() {
    if (this->fp.is_open()) this->fp.close();
}


EncBitstreamWriter &EncBitstreamWriter::operator+=(const EncBitstreamWriter &temp) {
    return dynamic_cast<EncBitstreamWriter &>(EncBitstreamVLC::operator+=(temp));
}

EncBitstreamWriter &EncBitstreamWriter::operator+=(const EncBitstreamVLC &temp) {
    return dynamic_cast<EncBitstreamWriter &>(EncBitstreamVLC::operator+=(temp));
}

uint EncBitstreamWriter::getTotalBytes() const {
    return totalBytes;
}
