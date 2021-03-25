#ifndef LIGHT_FIELD_CODEC_ENCSYMBOL_H
#define LIGHT_FIELD_CODEC_ENCSYMBOL_H

#include "EncBitstreamBuffer.h"

#include <math.h>
#include <vector>
#include <string>

typedef struct Symbol {
    int value;
    int len;
    unsigned int bitpattern;
}symbol;

class EncSymbol : public EncBitstreamBuffer{
public:
    EncSymbol(uint bufferSize = 50);

    int encodeLast(int last);
    int encodeRun(std::vector<int> run);
    int encodeRem(int rem);

    int encodeSymbol(std::vector<int> &e_buffer, int code, std::string type);

    int writeSyntaxElement(int info, int len);
    ~EncSymbol();

private:
    int expGolombEncode_ui(int value);
    void writeCode2Buffer(symbol *sym);
};


#endif //LIGHT_FIELD_CODEC_ENCSYMBOL_H
