//
// Created by igor on 10/04/2022.
//

#ifndef LF_CODEC_ENCODESYMBOLSMODEL_H
#define LF_CODEC_ENCODESYMBOLSMODEL_H

#include <cmath>
#include <vector>

#include "EncBitstreamBuffer.h"
#include "EntropyStructsModel.h"

using namespace std;

typedef unsigned long int CodeValue;                            // Type of arithmetic code value

#define Code_value_bits	32		                                     // Number of bits in a code value
#define Max_Value (((unsigned long long)1<<Code_value_bits)-1)   // Maximum value in a integer of 32 bits
#define First_qtr ((unsigned long int)Max_Value/4 + 1)           // Point after first quarter
#define Half ((unsigned long int)2*First_qtr)                    // Point after second quarter
#define Third_qtr ((unsigned long int)3*First_qtr)               // point after third quarter
#define Max_freq 131071

typedef struct Symbol {
    int value;
    int len;
    unsigned int bitpattern;
};

typedef struct {
    int number_of_symbols;
    int *index_to_symbol;
    int *symbol_to_index;
    unsigned long int *frequency;
    unsigned long int *cumulative_frequency;
}Encoder_Probabilistic_Model;

class EncodeSymbolsModel: public EncBitstreamBuffer{
public:
  EncodeSymbolsModel(uint bufferSize = 50);
    void encodeLREVector(std::vector<LRE_Struct> &lre);
    void encodeFinalBits();

    //ARITHMETIC ENCODER
    void arithEncodeSymbol(int symbol, int model);
    void arithDoneEncoding(void);
    void arithReset(void);

    //PROBABILISTIC MODELS
    int arithAddModel(void);
    void arithStartModel(int number_of_symbols, int model);
    void arithUpdateModel(int symbol, int m);
    void arithRestartModel(int m);
    void arithPrintModel(int model);
private:
    void writeCode2Buffer(Symbol *sym);
    void encodeExpGolomb(int symbol);
    void encodeBit(int bit);

    void arithBitPlusFollow(int bit);

    //CURRENT STATE OF THE ENCODING
    CodeValue low, high;                                       // End values in a range of current code-region
    long bits_to_follow;                                       // Number of opposite bits to output after the next bit.

    //PROBABILISTIC MODELS
    int number_of_models;
    Encoder_Probabilistic_Model *model;
};


#endif // LF_CODEC_ENCODESYMBOLSMODEL_H
