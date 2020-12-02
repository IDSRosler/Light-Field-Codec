#include "EncSymbol.h"

EncSymbol::EncSymbol(uint bufferSize) : EncBitstreamBuffer(bufferSize){}

EncSymbol::~EncSymbol() {
}

void EncSymbol::writeCode2Buffer(symbol *sym) {
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

int EncSymbol::encodeLast(int last) { // Max value 255 - 8 bits
    symbol sym{};

    sym.value = last;
    sym.len = 8;
    sym.bitpattern = last;

    this->writeCode2Buffer(&sym);

    return sym.len;
}

int EncSymbol::encodeSymbol(std::vector<int> & e_buffer, int code, std::string type){
    int m;
    symbol sym{};


    if (type == "8-bits"){ // 8 bits
        sym.value = code;
        sym.len = 8;
        sym.bitpattern = code;
    } else { // expgolomb code
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
            e_buffer.push_back(1);
        } else{
            e_buffer.push_back(0);
        }
        mask >>= 1u;
    }
}

int EncSymbol::encodeRun(std::vector<int> run) {
    int total_bits = 0;
    for (int i = 0; i < run.size(); ++i) {
        total_bits += this->expGolombEncode_ui(run[i]);
    }
    return total_bits;
}

int EncSymbol::encodeRem(int rem) {
    return this->expGolombEncode_ui(rem);
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