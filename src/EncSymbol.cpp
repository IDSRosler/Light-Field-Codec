#include "EncSymbol.h"

EncSymbol::EncSymbol(uint bufferSize) : EncBitstreamBuffer(bufferSize){}

EncSymbol::~EncSymbol() {
}

void EncSymbol::writeCode2Buffer(symbol *sym) {
    unsigned int mask = 1 << (sym->len - 1);
    Byte *byte_buf = &this->byte_buf;
    uint *bits_to_go = &this->bits_to_go;
    int i;

    if (sym->len < 33) {
        for (i = 0; i < sym->len; i++) {
            *byte_buf <<= 1u;

            if (sym->bitpattern & mask)
                *byte_buf |= 1u;

            mask >>= 1u;

            if ((--(*bits_to_go)) == 0) {
                *bits_to_go = 8;
                this->buffer[this->byte_pos++] = *byte_buf;
                *byte_buf = 0;
            }
        }
    } else { /*se->len >= 33*/
        // zeros
        for (i = 0; i < (sym->len - 32); i++) {
            *byte_buf <<= 1u;

            if ((--(*bits_to_go)) == 0) {
                *bits_to_go = 8;
                this->buffer[this->byte_pos++] = *byte_buf;
                *byte_buf = 0;
            }
        }
        // actual info
        mask = 1u << 31u;
        for (i = 0; i < 32; i++) {
            *byte_buf <<= 1u;

            if (sym->bitpattern & mask)
                *byte_buf |= 1u;

            mask >>= 1u;

            if ((--(*bits_to_go)) == 0) {
                *bits_to_go = 8;
                this->buffer[this->byte_pos++] = *byte_buf;
                *byte_buf = 0;
            }
        }
    }

    /*for (i = 0; i < sym->len; i++){
        *byte_buf <<= 1u;

        if (sym->bitparttern & mask){
            *byte_buf |= 1u;
        }

        mask >>= 1u;

        if ((--(*bits_to_go)) == 0) {
            *bits_to_go = 8;
            this->buffer[this->byte_pos++] = *byte_buf;
            *byte_buf = 0;
        }
    }*/
}

void EncSymbol::encodeLast(int last) { // Max value 255 - 8 bits
    symbol sym{};

    sym.value = last;
    sym.len = 8;
    sym.bitpattern = last;

    this->writeCode2Buffer(&sym);
}

void EncSymbol::encodeRun(std::vector<int> run) {
    for (int i = 0; i < run.size(); ++i) {
        this->expGolombEncode_ui(run[i]);
    }
}

int EncSymbol::expGolombEncode_ui(int value) {
    int m = 0, inf = 0;
    symbol sym{};

    sym.value = value;

    m = floor(log2(sym.value + 1));
    inf = sym.value + 1 - pow(2, m);

    sym.len = (m * 2) + 1;
    sym.bitpattern = sym.value + 1;

    this->writeCode2Buffer(&sym);

    return sym.len;
}

int EncSymbol::writeSyntaxElement(int info, int len) {
    symbol sym{};

    sym.value = info;
    sym.len = len;
    sym.bitpattern = info;

    this->writeCode2Buffer(&sym);

    return (sym.len);
}