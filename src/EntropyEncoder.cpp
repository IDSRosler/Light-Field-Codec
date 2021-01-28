#include "EntropyEncoder.h"

#include <utility>

EntropyEncoder::EntropyEncoder(EncoderParameters *parameters, uint bufferSize) : EncSymbol(bufferSize){
    this->parameters = parameters;
    this->open_file(this->parameters->getPathOutput() + "LightField.bin");

    this->report.openFiles(this->parameters->getPathOutput());
    this->report.setHeaders();

    /*//last update
    this->freqFile.open(this->parameters->getPathOutput() + "frequency.csv");
    this->freqFile << "LFBPU, Freq_sig_0, Freq_sig_1, Freq_grOne_0, Freq_grOne_1, Freq_grTwo_0, Freq_grTwo_1, Freq_sign_0, Freq_sign_1" << endl;*/
}

EntropyEncoder::~EntropyEncoder() = default;

void EntropyEncoder::encodeHypercube(int *bitstream, const Point4D &dim_block, int hypercube_pos, string channel) {
    this->hypercube = hypercube_pos;
    this->ch = std::move(channel);
    this->sig_sub = 0;
    this->n_sig_sub = 0;
    this->sig_coeff = 0;
    this->n_sig_coeff = 0;
    this->one = 0;
    this->two = 0;
    this->gr_two = 0;
    this->max_value = 0;
    this->mean_value = 0;

    int last_block = -1;
    vector<int> run;
    vector<SyntacticElements> lfbpu_elements; // light field base processing unit

    this->root = this->tree.CreateRoot(bitstream, dim_block);    // create tree
    this->tree.CreateTree(this->root, dim_block, {0,0,0,0});

    this->sig_coeff = this->root->att->sig_coeff;
    this->n_sig_coeff = static_cast<int>(dim_block.x * dim_block.y * dim_block.u * dim_block.v) - this->sig_coeff;
    this->one = this->root->att->n_one;
    this->two = this->root->att->n_two;
    this->gr_two = this->root->att->n_greater_than_two;
    this->max_value = this->root->att->max_value;
    this->mean_value = this->root->att->mean_value;

    this->model_8bits = this->arith_encoder.Add_model();
    this->model_exp = this->arith_encoder.Add_model();
    this->model_sig = this->arith_encoder.Add_model();
    this->model_grone = this->arith_encoder.Add_model();
    this->model_grtwo = this->arith_encoder.Add_model();
    this->model_sign = this->arith_encoder.Add_model();

    vector<int> freq;
    this->ComputeFrequency(freq);

    this->arith_encoder.Start_model(2, freq, this->model_8bits);
    this->arith_encoder.Start_model(2, freq, this->model_exp);
    this->arith_encoder.Start_model(2, freq, this->model_sig);
    this->arith_encoder.Start_model(2, freq, this->model_grone);
    this->arith_encoder.Start_model(2, freq, this->model_grtwo);
    this->arith_encoder.Start_model(2, freq, this->model_sign);

    this->report.setAtt(this->hypercube, this->ch);

    this->tree.ComputeLast(last_block);   // compute last (block level)

    this->tree.subpartitionReport(this->report);

    this->last = last_block;

    this->encodeSymbol(last_block, this->model_8bits, "8-bits");
    //this->encodeLast(last_block);

    this->tree.ComputeRun(run, this->last);  // compute run (block level)

    for (auto value : run){
        this->encodeSymbol(value, this->model_exp, "exp");
    }
    //this->encodeRun(run);

    run.clear();

    this->tree.ComputeSyntacticElements(lfbpu_elements, this->last);  // compute syntactic elements (coefficients level)

    this->sig_sub = lfbpu_elements.size();
    this->n_sig_sub = 256 - this->sig_sub;

    this->encodeSyntacticElements(lfbpu_elements);
    //this->EncodeSyntacticElements(lfbpu_elements);

    lfbpu_elements.clear();

    this->arith_encoder.Done_encoding();
    this->arith_encoder.Reset();
    this->e_buffer.clear();

    this->report.writeStatistics(this->last, this->sig_sub, this->n_sig_sub,
                                 this->sig_coeff, this->n_sig_coeff, this->one, this->two, this->gr_two,
                                 this->max_value, this->mean_value);

    this->tree.DeleteTree(&this->root); // delete tree
}

int EntropyEncoder::encodeSymbol(int code, int model, const std::string& type){
    symbol sym{};

    if (type == "bit"){ // bit
        sym.value = code;
        sym.len = 1;
        sym.bitpattern = code;
    } else if (type == "8-bits"){ // 8 bits
        sym.value = code;
        sym.len = 8;
        sym.bitpattern = code;
    } else { // expgolomb code
        int m;
        sym.value = code;
        m = floor(log2(sym.value + 1));
        sym.len = (m * 2) + 1;
        sym.bitpattern = sym.value + 1;
    }

    // write buffer
    unsigned int mask = 1 << (sym.len - 1);
    int i;
    for (i = 0; i < sym.len; i++){
        if (sym.bitpattern & mask){
            this->arith_encoder.Encode_symbol(1, model);
            this->arith_encoder.update_model(1, model);
        } else{
            this->arith_encoder.Encode_symbol(0, model);
            this->arith_encoder.update_model(0, model);
        }
        mask >>= 1u;
    }

    return sym.len;
}

void EntropyEncoder::encodeSyntacticElements(vector<SyntacticElements> &lfbpu){
    // Add bits on buffer
    for (auto & i : lfbpu) {
        if (i.last >= 0){
            this->encodeSymbol(i.last, this->model_8bits, "8-bits"); // encode last
            for (int sig: i.sig){ // encode sig
                this->encodeSymbol(sig, this->model_sig, "bit");
            }
            for (int gr_one : i.gr_one){ // encode gr_one
                this->encodeSymbol(gr_one, this->model_grone, "bit");
            }
            for (int gr_t : i.gr_two){ // encode gr_two
                this->encodeSymbol(gr_t, this->model_grtwo, "bit");
            }
            for (int sign : i.sign){ // encode sign
                this->encodeSymbol(sign, this->model_sign, "bit");
            }
            for (int rem : i.rem){ // encode rem
                this->encodeSymbol(rem, this->model_exp, "exp");
            }
        }
    }
}

void EntropyEncoder::EncodeSyntacticElements(vector<SyntacticElements> &lfbpu) {
    int sig_model = this->arith_encoder.Add_model();
    int gr_one_model = this->arith_encoder.Add_model();
    int gr_two_model = this->arith_encoder.Add_model();
    int sign_model = this->arith_encoder.Add_model();

    ElementsFrequency elem_freq;

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

    for (auto & i : lfbpu) {
        if (i.last > -1){
            this->encodeLast(i.last); // encode last

            for (auto sig: i.sig){ // encode sig
                this->arith_encoder.Encode_symbol(sig, sig_model);
            }
            for (auto gr_one : i.gr_one){ // encode gr_one
                this->arith_encoder.Encode_symbol(gr_one, gr_one_model);
            }
            for (auto gr_t : i.gr_two){ // encode gr_two
                this->arith_encoder.Encode_symbol(gr_t, gr_two_model);
            }
            for (auto sign : i.sign){ // encode sign
                this->arith_encoder.Encode_symbol(sign, sign_model);
            }
            this->arith_encoder.Done_encoding();

            for (auto rem : i.rem){ // encode rem
                this->encodeRem(rem);
            }
        }
    }
    this->arith_encoder.Reset();
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

    this->report.closeFiles();

   /* //last update
    if (this->freqFile.is_open()) this->freqFile.close();*/
}

uint EntropyEncoder::getTotalBytes() const {
    return totalBytes;
}

void EntropyEncoder::open_file(const string &filename) {
    this->outputFile.open(filename, std::ios::binary);
    assert(this->outputFile.is_open());
}

void EntropyEncoder::ComputeFrequency(vector<int> &freq){
    int qt_0 = 0, f_0;

    for (int i : this->e_buffer) {
        if (i == 0){
            ++qt_0;
        }
    }

    f_0 = std::round(( (float)qt_0 / this->e_buffer.size() ) * 100);

    freq.push_back(f_0);
    freq.push_back(100 - f_0);
}

void EntropyEncoder::ComputeFrequency(vector<SyntacticElements> &lfbpu, ElementsFrequency& elem_freq) {
    int size, q0 = 0,
        cont_sig = 0, cont_one = 0,
        cont_two = 0,cont_sign = 0,
        sig0 = 0, sum_sig0 = 0,
        one0 = 0, sum_one0 = 0,
        two0 = 0, sum_two0 = 0,
        sign0 = 0, sum_sing0 = 0;

    for (auto & i : lfbpu) {

        size = i.sig.size();
        if (size > 0){
            for (int j = 0; j < size; ++j) {    //freq sig
                if (i.sig[j] == 0)
                    ++q0;
            }
            sig0 = std::round(((float)q0/size)*100);
            sum_sig0 += sig0;
            ++cont_sig;
            q0 = 0;
        }

        size = i.gr_one.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq one
                if (i.gr_one[j] == 0)
                    ++q0;
            }
            one0 = std::round(((float)q0/size)*100);
            sum_one0 += one0;
            ++cont_one;
            q0 = 0;
        }

        size = i.gr_two.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq tow
                if (i.gr_two[j] == 0)
                    ++q0;
            }
            two0 = std::round(((float)q0/size)*100);
            sum_two0 += two0;
            ++cont_two;
            q0 = 0;
        }

        size = i.sign.size();
        if (size > 0) {
            for (int j = 0; j < size; ++j) {    //freq sign
                if (i.sign[j] == 0)
                    ++q0;
            }
            sign0 = std::round(((float)q0/size)*100);
            sum_sing0 += sign0;
            ++cont_sign;
            q0 = 0;
        }

       /* this->freqFile << i << "," <<
                           sig0 << "," <<
                           100 - sig0 << "," <<
                           one0 << "," <<
                           100 - one0 << "," <<
                           two0 << "," <<
                           100 - two0 << "," <<
                           sign0 << "," <<
                           100 - sign0 << "," << endl;
*/
    }

    int sig_mean_0 = std::round(sum_sig0/cont_sig);
    int sig_mean_1 = 100 - sig_mean_0;
    int one_mean_0 = std::round(sum_one0/cont_one);
    int one_mean_1 = 100 - one_mean_0;
    int two_mean_0 = std::round(sum_two0/cont_two);
    int two_mean_1 = 100 - two_mean_0;
    int sign_mean_0 = std::round(sum_sing0/cont_sign);
    int sign_mean_1 = 100 - sign_mean_0;

    elem_freq.setFrequency(sig_mean_0, sig_mean_1, one_mean_0, one_mean_1, two_mean_0, two_mean_1, sign_mean_0, sign_mean_1);

    /*this->freqFile << endl;
    this->freqFile << "Mean_of_Frequency_per_Hypercube, " <<
        sig_mean_0 << "," <<
        sig_mean_1 << "," <<
        one_mean_0 << "," <<
        one_mean_1 << "," <<
        two_mean_0 << "," <<
        two_mean_1 << "," <<
        sign_mean_0 << "," <<
        sign_mean_1 << endl;
    this->freqFile << endl;*/
}
