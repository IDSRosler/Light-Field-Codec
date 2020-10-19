#include "EntropyEncoder.h"

EntropyEncoder::EntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncSymbol(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");

    //last update
    this->freqFile.open(this->parameters->getPathOutput() + "frequency.csv");
    this->freqFile << "LFBPU, Freq_sig_0, Freq_sig_1, Freq_grOne_0, Freq_grOne_1, Freq_grTwo_0, Freq_grTwo_1, Freq_sign_0, Freq_sign_1" << endl;
}

EntropyEncoder::~EntropyEncoder() {
}

void EntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block) {
    int last = -1;
    vector<int> run;
    vector<SyntacticElements> lfbpu_elements; // light field base processing unit

    this->root = this->tree.CreateRoot(bitstream, dim_block);    // create tree
    this->tree.CreateTree(this->root, dim_block, {0,0,0,0});

    this->tree.ComputeLast(last);   // compute last (block level)
    this->encodeLast(last);

    if (last > 0){
        this->tree.ComputeRun(run, last);  // compute run (block level)
        this->encodeRun(run);
        run.clear();

        this->tree.ComputeSyntacticElements(lfbpu_elements, last);  // compute syntactic elements (coefficients level)

        // last update
        this->ComputeFrequency(lfbpu_elements);

        this->EncodeSyntacticElements(lfbpu_elements);

        lfbpu_elements.clear();
    }

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
        if (lfbpu[i].last != -1){
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
            for (auto rem : lfbpu[i].rem){ // encode rem
                this->encodeRem(rem);
            }
        }
    }
}

void EntropyEncoder::write_completedBytes() {
    if (this->byte_pos == 0) return;
    this->outputFile.write((char *) this->buffer, this->byte_pos);
    this->totalBytes += this->byte_pos;
    this->byte_pos = 0;
}

void EntropyEncoder::finish_and_write() {
    if(this->bits_to_go < 8)
        this->writeSyntaxElement(0, this->bits_to_go);
    this->write_completedBytes();
    if (this->outputFile.is_open()) this->outputFile.close();

    //last update
    if (this->freqFile.is_open()) this->freqFile.close();
}

uint EntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

void EntropyEncoder::open_file(const string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}

void EntropyEncoder::ComputeFrequency(vector<SyntacticElements> lfbpu) {
    int size = 0,
        q0 = 0,
        sig0 = 0,
        sig1 = 0,
        one0 = 0,
        one1 = 0,
        two0 = 0,
        two1 = 0,
        sign0 = 0,
        sign1 = 0;

    for (int i = 0; i < lfbpu.size(); ++i) {

        size = lfbpu[i].sig.size();
        if (size > 0){
            for (int j = 0; j < size; ++j) {    //freq sig
                if (lfbpu[i].sig[j] == 0)
                    ++q0;
            }
            sig0 = round((float(q0)/size)*100);
            sig1 = 100 - sig0;
            q0 = 0;
        }

        size = lfbpu[i].gr_one.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq one
                if (lfbpu[i].gr_one[j] == 0)
                    ++q0;
            }
            one0 = round((float(q0)/ size) * 100);
            one1 = 100 - one0;
            q0 = 0;
        }

        size = lfbpu[i].gr_two.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq tow
                if (lfbpu[i].gr_two[j] == 0)
                    ++q0;
            }
            two0 = round((float(q0)/ size) * 100);
            two1 = 100 - two0;
            q0 = 0;
        }

        size = lfbpu[i].sign.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq sign
                if (lfbpu[i].sign[j] == 0)
                    ++q0;
            }
            sign0 = round((float(q0)/ size) * 100);
            sign1 = 100 - sign0;
            q0 = 0;
        }

        this->freqFile << i << "," <<
                           sig0 << "," <<
                           sig1 << "," <<
                           one0 << "," <<
                           one1 << "," <<
                           two0 << "," <<
                           two1 << "," <<
                           sign0 << "," <<
                           sign1 << "," << endl;

    }
}