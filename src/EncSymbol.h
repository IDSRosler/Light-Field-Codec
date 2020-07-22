#ifndef LIGHT_FIELD_CODEC_ENCSYMBOL_H
#define LIGHT_FIELD_CODEC_ENCSYMBOL_H

#include "EncBitstreamBuffer.h"

#include <math.h>
#include <vector>

typedef struct Symbol {
    int value;
    int len;
    unsigned int bitpattern;
}symbol;

class EncSymbol : public EncBitstreamBuffer{
public:
    EncSymbol(uint bufferSize = 50);

    void encodeLast(int last);
    void encodeRun(std::vector<int> run);

    int writeSyntaxElement(int info, int len);
    ~EncSymbol();

private:
    int expGolombEncode_ui(int value);
    void writeCode2Buffer(symbol *sym);
};


#endif //LIGHT_FIELD_CODEC_ENCSYMBOL_H
