//
// Created by igor on 10/04/2022.
//

#ifndef LF_CODEC_ENCODESYMBOL_H
#define LF_CODEC_ENCODESYMBOL_H

#include <cmath>
#include <vector>

#include "EncBitstreamBuffer.h"
#include "ArithmeticStructures.h"

typedef struct Symbol {
    int value;
    int len;
    unsigned int bitpattern;
};


class EncodeSymbol : public EncBitstreamBuffer{
public:
    EncodeSymbol(uint bufferSize = 50);
    void encodeLREVector(std::vector<LRE_Struct> &lre);
    void encodeFinalBits();
private:
    void writeCode2Buffer(Symbol *sym);
    void encodeExpGolomb(int symbol);
    void encodeBit(int bit);
};


#endif //LF_CODEC_ENCODESYMBOL_H
