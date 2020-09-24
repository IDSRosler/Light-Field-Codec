#ifndef LIGHT_FIELD_CODEC_ARITHMETICENCODER_H
#define LIGHT_FIELD_CODEC_ARITHMETICENCODER_H

#include "EncSymbol.h"

#include <fstream>
#include <iostream>

using namespace std;

typedef unsigned long int Code_Value;                           // Type of a arithmetic code value
/*typedef unsigned char Byte;                                     // Type byte to buffer*/

#define Code_value_bits	32		                                // Number of bits in a code value
#define Max_Value (((unsigned long long)1<<Code_value_bits)-1)  // Maximum value in a integer of 32 bits
#define First_qtr ((unsigned long int)Max_Value/4 + 1)       // Point after first quarter
#define Half ((unsigned long int)2*First_qtr)                   // Point after second quarter
#define Third_qtr ((unsigned long int)3*First_qtr)              // point after thrird quarter

typedef struct {
    int number_of_symbols;
    int *index_to_symbol;
    int *symbol_to_index;
    unsigned long int *frequency;
    unsigned long int *cumulative_frequency;
}Encoder_Probabilistic_Model;

class ArithmeticEncoder : EncSymbol{
public:
    ArithmeticEncoder(void);

    //ARITHMETIC ENCODER
    void Encode_symbol(int symbol, int model);
    void Done_encoding(void);
    void Set_output_file(char *path);
    void Reset(void);

    //PROBABILISTIC MODELS
    int Add_model(void);
    void Start_model(int number_of_symbols, vector<int> frequency,int model);
    void Print_model(int model);

private:
    void Bit_plus_follow(int bit);
    void Output_bit(int bit);
    void Done_output_bits(void);

    /*//OUTPUT FILE
    ofstream file;*/

    //CURRENT STATE OF THE ENCODING
    Code_Value low, high;                                       // Ends values in a range of current code-region
    long bits_to_follow;                                        // Number of opposite bits to output after the next bit.

    /*//THE BIT BUFFER
    Byte buffer;                                                // Bits buffered for output
    int bits_to_go;                                             // Number of bits still in buffer*/

    //PROBABILISTIC MODELS
    int number_of_models;
    Encoder_Probabilistic_Model *model;
};

#endif //LIGHT_FIELD_CODEC_ARITHMETICENCODER_H
