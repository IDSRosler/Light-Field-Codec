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

    //Todo: Fix last read and the buffer write

    this->tree.ComputeRun(run, last);  // compute run (block level)
    this->encodeRun(run);
    run.clear();

    //this->write_completedBytes();

    this->tree.ComputeSyntacticElements(lfbpu_elements, last);  // compute syntactic elements (coefficients level)

    this->EncodeSyntacticElements(lfbpu_elements);

    lfbpu_elements.clear();
    this->tree.DeleteTree(&this->root); // delete tree
}

void EntropyEncoder::EncodeSyntacticElements(vector<SyntacticElements> lfbpu) {
    int sig_model = this->arith_encoder.Add_model();
    int gr_one_model = this->arith_encoder.Add_model();
    int gr_two_model = this->arith_encoder.Add_model();
    int sign_model = this->arith_encoder.Add_model();

    ElementsFrequency elem_freq;
    elem_freq.setFrequency();

    this->arith_encoder.Start_model(2, elem_freq.sig, sig_model);
    this->arith_encoder.Start_model(2, elem_freq.gr_one, gr_one_model);
    this->arith_encoder.Start_model(2, elem_freq.gr_two, gr_two_model);
    this->arith_encoder.Start_model(2, elem_freq.sign, sign_model);

    /*this->arith_encoder.Print_model(sig_model);
    this->arith_encoder.Print_model(gr_one_model);
    this->arith_encoder.Print_model(gr_two_model);
    this->arith_encoder.Print_model(sign_model);*/

    for (int i = 0; i < lfbpu.size(); ++i) {
        this->encodeLast(lfbpu[i].last); // encode last
        for (auto sig: lfbpu[i].sig){ // encode sig
            this->arith_encoder.Encode_symbol(sig, sig_model);
        }
        for (auto gr_one : lfbpu[i].gr_one){ // encode gr_one
            this->arith_encoder.Encode_symbol(gr_one, gr_one_model);
        }
        for (auto gr_two : lfbpu[i].gr_two){ // encode gr_two
            this->arith_encoder.Encode_symbol(gr_two, gr_two_model);
        }
        for (auto sign : lfbpu[i].sign){ // encode sign
            this->arith_encoder.Encode_symbol(sign, sign_model);
        }
        this->arith_encoder.Done_encoding();

        //Todo: encode rem here
    }
}

void EntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->writeSyntaxElement(0, this->bits_to_go);
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