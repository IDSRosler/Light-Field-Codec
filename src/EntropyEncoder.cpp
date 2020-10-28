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

    // last update
    this->ComputeFrequency(lfbpu, elem_freq);

    //elem_freq.setFrequency(50,50,50,50,50,50,50,50);

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

void EntropyEncoder::ComputeFrequency(vector<SyntacticElements> lfbpu, ElementsFrequency& elem_freq) {
    int size, q0 = 0,
        cont_sig = 0, cont_one = 0,
        cont_two = 0,cont_sign = 0,
        sig0 = 0, sum_sig0 = 0,
        one0 = 0, sum_one0 = 0,
        two0 = 0, sum_two0 = 0,
        sign0 = 0, sum_sing0 = 0;

    for (int i = 0; i < lfbpu.size(); ++i) {

        size = lfbpu[i].sig.size();
        if (size > 0){
            for (int j = 0; j < size; ++j) {    //freq sig
                if (lfbpu[i].sig[j] == 0)
                    ++q0;
            }
            sig0 = round((float(q0)/size)*100);
            sum_sig0 += sig0;
            ++cont_sig;
            q0 = 0;
        }

        size = lfbpu[i].gr_one.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq one
                if (lfbpu[i].gr_one[j] == 0)
                    ++q0;
            }
            one0 = round((float(q0)/ size) * 100);
            sum_one0 += one0;
            ++cont_one;
            q0 = 0;
        }

        size = lfbpu[i].gr_two.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq tow
                if (lfbpu[i].gr_two[j] == 0)
                    ++q0;
            }
            two0 = round((float(q0)/ size) * 100);
            sum_two0 += two0;
            ++cont_two;
            q0 = 0;
        }

        size = lfbpu[i].sign.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq sign
                if (lfbpu[i].sign[j] == 0)
                    ++q0;
            }
            sign0 = short(round((float(q0)/ size) * 100));
            sum_sing0 += sign0;
            ++cont_sign;
            q0 = 0;
        }

        this->freqFile << i << "," <<
                           sig0 << "," <<
                           100 - sig0 << "," <<
                           one0 << "," <<
                           100 - one0 << "," <<
                           two0 << "," <<
                           100 - two0 << "," <<
                           sign0 << "," <<
                           100 - sign0 << "," << endl;

    }

    int sig_mean_0 = round(float(sum_sig0)/cont_sig);
    int sig_mean_1 = 100 - sig_mean_0;
    int one_mean_0 = round(float(sum_one0)/cont_one);
    int one_mean_1 = 100 - one_mean_0;
    int two_mean_0 = round(float(sum_two0)/cont_two);
    int two_mean_1 = 100 - two_mean_0;
    int sign_mean_0 = round(float(sum_sing0)/cont_sign);
    int sign_mean_1 = 100 - sign_mean_0;

    elem_freq.setFrequency(sig_mean_0, sig_mean_1, one_mean_0, one_mean_1, two_mean_0, two_mean_1, sign_mean_0, sign_mean_1);

    this->freqFile << endl;
    this->freqFile << "Mean_of_Frequency_per_Hypercube, " <<
        sig_mean_0 << "," <<
        sig_mean_1 << "," <<
        one_mean_0 << "," <<
        one_mean_1 << "," <<
        two_mean_0 << "," <<
        two_mean_1 << "," <<
        sign_mean_0 << "," <<
        sign_mean_1 << endl;
    this->freqFile << endl;
}