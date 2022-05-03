//
// Created by igor on 10/04/2022.
//

#include "EncodeSymbolsModel.h"

EncodeSymbolsModel::EncodeSymbolsModel(uint bufferSize) : EncBitstreamBuffer(bufferSize){
    this->low = 0;
    this->high = Max_Value;
    this->bits_to_follow = 0;

    this->number_of_models = 0;
    this->model = nullptr;
}

// ----------------------------------- LRE Encode --------------------------------------------------
void EncodeSymbolsModel::encodeLREVector(std::vector<LRE_Struct> &lre) {
    int size = lre.size();
    if (size >= 0) {
        this->encodeExpGolomb(size);
    }
    while (!lre.empty()){
        LRE_Struct lreValue = lre.back();
        this->encodeBit(lreValue.level);
        this->encodeExpGolomb(lreValue.run);
        lre.pop_back();
    }
}

// ----------------------------------- Exp-Golomb Encode -------------------------------------------
void EncodeSymbolsModel::encodeExpGolomb(int symbol) {
    Symbol sym{};
    int m;
    sym.value = symbol;
    m = floor(log2(sym.value + 1));
    sym.len = (m * 2) + 1;
    sym.bitpattern = sym.value + 1;

    writeCode2Buffer(&sym);
}

// --------------------------------------- Bits Level ----------------------------------------------
void EncodeSymbolsModel::encodeBit(int bit) {
    Symbol sym{};

    sym.value = bit;
    sym.len = 1;
    sym.bitpattern = bit;

    writeCode2Buffer(&sym);
}

void EncodeSymbolsModel::encodeFinalBits() {
    for (int i = 0; i < this->bits_to_go; i++){
        encodeBit(0);
    }
    this->bits_to_go = 0;
}

// --------------------------------- Write int the Buffer ------------------------------------------
void EncodeSymbolsModel::writeCode2Buffer(Symbol *sym) {
    unsigned int mask = 1 << (sym->len - 1);
    Byte *byte_buf = &this->byte_buf;
    uint *bits_to_go = &this->bits_to_go;
    int i;

    for (i = 0; i < sym->len; i++){
        *byte_buf <<= 1u;

        if (sym->bitpattern & mask){
            *byte_buf |= 1u;
        }

        mask >>= 1u;

        if ((--(*bits_to_go)) == 0) {
            *bits_to_go = 8;
            this->buffer[this->byte_pos++] = *byte_buf;
            *byte_buf = 0;
        }
    }
}

// ------------------------------------- Arithmetic Encode -----------------------------------------
void EncodeSymbolsModel::arithEncodeSymbol(int symbol, int model) {
    long long range;        // Size of the current_code region

    symbol = this->model[model].symbol_to_index[symbol];

    range = (long long)(this->high - this->low)+1;

    this->high = this->low +
                 (range * this->model[model].cumulative_frequency[symbol-1]) /
                 this->model[model].cumulative_frequency[0] - 1;

    this->low = this->low +
                (range * this->model[model].cumulative_frequency[symbol]) /
                this->model[model].cumulative_frequency[0];

    while (true){
        if (this->high < Half){
            this->arithBitPlusFollow(0);
        } else if (this->low >= Half){
            this->arithBitPlusFollow(1);
            this->low = this->low - Half;
            this->high = this->high - Half;
        } else if (this->low >= First_qtr && this->high < Third_qtr){
            ++this->bits_to_follow;
            this->low = this->low - First_qtr;
            this->high = this->high - First_qtr;
        } else break;
        this->low = 2 * this->low;
        this->high = 2 * this->high + 1;
    }
}

void EncodeSymbolsModel::arithDoneEncoding() {
    ++this->bits_to_follow;
    if (this->low < First_qtr) this->arithBitPlusFollow(0);
    else this->arithBitPlusFollow(1);

    this->encodeFinalBits();
}

void EncodeSymbolsModel::arithBitPlusFollow(int bit) {
    this->encodeBit(bit);
    while (this->bits_to_follow > 0){
        this->encodeBit(!bit);
        --this->bits_to_follow;
    }
}

void EncodeSymbolsModel::arithReset() {
    this->low = 0;
    this->high = Max_Value;
    this->bits_to_follow = 0;

    this->buffer = 0;
    this->bits_to_go = 8;

    for (int i = 0; i < this->number_of_models; ++i) {
        if (this->model[i].cumulative_frequency != nullptr) free(this->model[i].cumulative_frequency);
        if (this->model[i].frequency!= nullptr) free(this->model[i].frequency);
        if (this->model[i].symbol_to_index != nullptr) free(this->model[i].symbol_to_index);
        if (this->model[i].index_to_symbol != nullptr) free(this->model[i].index_to_symbol);
    }

    this->number_of_models = 0;
    if (this->model != nullptr) free(this->model);
    this->model = nullptr;
}

/*
****************************************************************************************************************
                                            ADAPTIVE MODEL
 ****************************************************************************************************************
 */

int EncodeSymbolsModel::arithAddModel() {
    ++this->number_of_models;
    if (this->model == nullptr) {
        if ((this->model = (Encoder_Probabilistic_Model *) malloc(this->number_of_models*sizeof(Encoder_Probabilistic_Model))) ==
            nullptr){ cout << "Memory allocation error" << endl; exit(1);}
    } else {
        if ((this->model = (Encoder_Probabilistic_Model *) realloc(this->model, this->number_of_models*sizeof(Encoder_Probabilistic_Model))) ==
            nullptr) { cout << "Memory allocation error" << endl; exit(1); }
    }
    this->model[this->number_of_models-1].number_of_symbols = 0;

    return this->number_of_models - 1;
}

void EncodeSymbolsModel::arithStartModel(int number_of_symbols, int model) {
    this->model[model].number_of_symbols = number_of_symbols;

    // CUMULATIVE FREQUENCY
    if ((this->model[model].cumulative_frequency =
         (unsigned long int *) malloc((number_of_symbols+1)*sizeof(unsigned long int))) ==
        nullptr) { cout << "Memory allocation error" << endl; exit(1); }

    // FREQUENCY
    if ((this->model[model].frequency =
         (unsigned long int *) malloc((number_of_symbols+1)*sizeof(unsigned long int))) ==
        nullptr) { cout << "Memory allocation error" << endl; exit(1); }

    // SYMBOL TO INDEX
    if ((this->model[model].symbol_to_index =
         (int *) malloc((number_of_symbols+1)*sizeof(int))) ==
        nullptr) { cout << "Memory allocation error" << endl; exit(1); }

    // INDEX TO SYMBOL
    if ((this->model[model].index_to_symbol =
         (int *) malloc((number_of_symbols+1)*sizeof(int))) ==
        nullptr) { cout << "Memory allocation error" << endl; exit(1); }

    this->model[model].index_to_symbol[0] = 0;

    for (int i = 0; i < number_of_symbols; ++i) {
        this->model[model].symbol_to_index[i] = i+1;
        this->model[model].index_to_symbol[i+1] = i;
    }

    this->model[model].cumulative_frequency[number_of_symbols] = 0;
    this->model[model].frequency[0] = 0;

    for (int j = number_of_symbols; j > 0; --j) {
        this->model[model].frequency[j] = 1;
        this->model[model].cumulative_frequency[j-1] = this->model[model].cumulative_frequency[j] +
                                                         this->model[model].frequency[j];
    }
}

void EncodeSymbolsModel::arithUpdateModel(int symbol, int m) {
    int i;
    unsigned long int cum;

    symbol = this->model[m].symbol_to_index[symbol];

    if (this->model[m].cumulative_frequency[0]==Max_freq)
    {
        cum = 0;

        for (i = this->model[m].number_of_symbols; i>=0; i--)
        {
            this->model[m].frequency[i] = (model[m].frequency[i]+1)/2;
            this->model[m].cumulative_frequency[i] = cum;
            cum += this->model[m].frequency[i];
        }
    }
    this->model[m].frequency[0] = 0;

    for (i = symbol; this->model[m].frequency[i]==this->model[m].frequency[i-1]; i--) ;
    if (i<symbol)
    {
        int ch_i, ch_symbol;

        ch_i = this->model[m].index_to_symbol[i];
        ch_symbol = this->model[m].index_to_symbol[symbol];
        this->model[m].index_to_symbol[i] = ch_symbol;
        this->model[m].index_to_symbol[symbol] = ch_i;
        this->model[m].symbol_to_index[ch_i] = symbol;
        this->model[m].symbol_to_index[ch_symbol] = i;
    }

    this->model[m].frequency[i] += 1;

    while (i>0)
    {
        i--;
        this->model[m].cumulative_frequency[i]++;
    }
}

void EncodeSymbolsModel::arithRestartModel(int m) {
    int number_of_symbols = this->model[m].number_of_symbols;
    for (int i = 0; i < number_of_symbols; ++i) {
        this->model[m].symbol_to_index[i] = i+1;
        this->model[m].index_to_symbol[i+1] = i;
    }

    this->model[m].cumulative_frequency[number_of_symbols] = 0;
    this->model[m].frequency[0] = 0;

    for (int j = number_of_symbols; j > 0; --j) {
        this->model[m].frequency[j] = 1;
        this->model[m].cumulative_frequency[j-1] = this->model[m].cumulative_frequency[j] +
                                                         this->model[m].frequency[j];
    }
}

void EncodeSymbolsModel::arithPrintModel(int model) {
    cout << "Frequency[model = " << model << "] = [";
    for (int i = 0; i <= this->model[model].number_of_symbols; ++i) {
        cout << " " << this->model[model].frequency[i];
    }
    cout << " ]" << endl;

    cout << "Cumulative Frequency[model = " << model << "] = [";
    for (int i = 0; i <= this->model[model].number_of_symbols; ++i) {
        cout << " " << this->model[model].cumulative_frequency[i];
    }
    cout << " ]" << endl;

    cout << "Symbol to index[model = " << model << "] = [";
    for (int i = 0; i <= this->model[model].number_of_symbols; ++i) {
        cout << " " << this->model[model].symbol_to_index[i];
    }
    cout << " ]" << endl;

    cout << "Index to symbol[model = " << model << "] = [";
    for (int i = 0; i <= this->model[model].number_of_symbols; ++i) {
        cout << " " << this->model[model].index_to_symbol[i];
    }
    cout << " ]\n" << endl;
}