//
// Created by igor on 10/04/2022.
//

#include "EncodeSymbol.h"

EncodeSymbol::EncodeSymbol(uint bufferSize) : EncBitstreamBuffer(bufferSize){}

void EncodeSymbol::encodeLREVector(std::vector<LRE_Struct> &lre) {
    while (!lre.empty()){
        LRE_Struct lreValue = lre.back();
        this->encodeBit(lreValue.level);
        this->encodeExpGolomb(lreValue.run);
        lre.pop_back();
    }
}

void EncodeSymbol::writeCode2Buffer(Symbol *sym) {
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

void EncodeSymbol::encodeExpGolomb(int symbol) {
    Symbol sym{};
    int m;
    sym.value = symbol;
    m = floor(log2(sym.value + 1));
    sym.len = (m * 2) + 1;
    sym.bitpattern = sym.value + 1;

    writeCode2Buffer(&sym);
}

void EncodeSymbol::encodeBit(int bit) {
    Symbol sym{};

    sym.value = bit;
    sym.len = 1;
    sym.bitpattern = bit;

    writeCode2Buffer(&sym);
}

void EncodeSymbol::encodeFinalBits() {
    for (int i = 0; i < this->bits_to_go; i++){
        encodeBit(0);
    }
    this->bits_to_go = 0;
}