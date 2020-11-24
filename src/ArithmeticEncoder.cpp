#include "ArithmeticEncoder.h"

/*
****************************************************************************************************************
                                        ARITHMETIC ENCODING
 ****************************************************************************************************************
 */

ArithmeticEncoder :: ArithmeticEncoder(Byte &buffer, uint &bits_to_go, Byte &byte_buf, uint &byte_pos) {
    this->low = 0;
    this->high = Max_Value;
    this->bits_to_follow = 0;

    this->buffer = &buffer;
    this->byte_pos = &byte_pos;
    this->bits_to_go = &bits_to_go;
    this->byte_buf = &byte_buf;

    this->local_buffer = 0;
    this->local_bits_to_go = 8;

    this->number_of_models = 0;
    this->model = nullptr;

    this->total_bits = 0;
}

//ENCODE A SYMBOL
void ArithmeticEncoder :: Encode_symbol(int symbol, int model) {
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
            this->Bit_plus_follow(0);
        } else if (this->low >= Half){
            this->Bit_plus_follow(1);
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

// END OF ENCODING
int ArithmeticEncoder::Done_encoding() { // Output two bits that select the quarter that the current code range contains
    int bits = this->total_bits;

    ++this->bits_to_follow;
    if (this->low < First_qtr) this->Bit_plus_follow(0);
    else this->Bit_plus_follow(1);

    this->Done_output_bits();

    this->total_bits = 0;

    return bits;
}

// OUTPUT BITS PLUS THE FOLLOWING OPPOSITE BITS
void ArithmeticEncoder :: Bit_plus_follow(int bit) {
    this->Output_bit(bit);
    while (this->bits_to_follow > 0){
        this->Output_bit(!bit);
        --this->bits_to_follow;
    }
}

// WRITE A BIT IN THE BUFFER
void ArithmeticEncoder :: Output_bit(int bit) {
    this->local_buffer >>= 1;               // Shift a bit to right
    if (bit) this->local_buffer |= 0x80;    // Put bit in top of buffer
    --this->local_bits_to_go;
    if (this->local_bits_to_go == 0){
        this->local_bits_to_go = 8;
        this->writeCode2Buffer();
    }
}

// WRITE THE LAST BITS
void ArithmeticEncoder :: Done_output_bits() {
    this->local_buffer = this->local_buffer >> this->local_bits_to_go;
    this->writeCode2Buffer();
}

void ArithmeticEncoder::writeCode2Buffer() {
    unsigned int mask = 1 << 7;
    int i;

    for (i = 0; i < 8; i++){
        *this->byte_buf <<= 1u;

        if (this->local_buffer & mask){
            *this->byte_buf |= 1u;
        }

        mask >>= 1u;

        if ((--(*this->bits_to_go)) == 0) {
            *this->bits_to_go = 8;
            this->buffer[(*this->byte_pos)++] = *this->byte_buf;
            *this->byte_buf = 0;
        }
    }
    this->local_bits_to_go = 8;
    this->local_buffer = 0;
    this->total_bits += 8;
}

// CODEC RESET
void ArithmeticEncoder :: Reset() {
    this->low = 0;
    this->high = Max_Value;
    this->bits_to_follow = 0;

    this->local_buffer = 0;
    this->local_bits_to_go = 8;

    this->total_bits = 0;

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

// ADD A NEW MODEL
int ArithmeticEncoder::Add_model() {
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

// INITIALIZE A MODEL
void ArithmeticEncoder::Start_model(int number_of_symbols, vector<int> frequency,int model) {
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

    for (int i = 0; i < frequency.size(); ++i) {
        this->model[model].frequency[i+1] = frequency[i];
    }

    for (int j = number_of_symbols; j > 0; --j) {
        //this->model[model].frequency[j] = 1;
        this->model[model].cumulative_frequency[j-1] = this->model[model].cumulative_frequency[j] +
                                                 this->model[model].frequency[j];
    }
}

void ArithmeticEncoder::Print_model(int model) {
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