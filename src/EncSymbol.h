#ifndef LIGHT_FIELD_CODEC_ENCSYMBOL_H
#define LIGHT_FIELD_CODEC_ENCSYMBOL_H

#include "EncBitstreamBuffer.h"
#include <math.h>

struct Symbol {
    int value;
    int len;
    unsigned int bitparttern;
};

class EncSymbol : public EncBitstreamBuffer{
public:
    EncSymbol(uint bufferSize = 50);
    ~EncSymbol();
private:
    int expGolombEncode_ui(int value);
    void writeCode2Buffer(Symbol *sym);
};


#endif //LIGHT_FIELD_CODEC_ENCSYMBOL_H
