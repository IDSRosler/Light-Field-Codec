#ifndef LF_CODEC_ENCBITSTREAMVLC_H
#define LF_CODEC_ENCBITSTREAMVLC_H

#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include <climits>
#include <cassert>
#include <cmath>
#include <bitset>


#include "Point4D.h"
#include "LRE.h"
#include "CAVLC.h"
#include "EncBitstreamBuffer.h"
#include "EncoderParameters.h"

#define TOTRUN_NUM       15
#define RUNBEFORE_NUM_M1  6

typedef struct Syntaxelement_enc {
    int type; //!< type of syntax element for data part.
    int value1; //!< numerical value of syntax element
    int value2; //!< for blocked symbols, e.g. run/level
    int len; //!< length of code
    int inf; //!< info part of UVLC code
    unsigned int bitpattern; //!< UVLC bitpattern

} syntaxelement_enc;

class EncBitstreamVLC : public EncBitstreamBuffer {

    LRE *lre;

    uint size, entropy_type{1};

public:

    EncBitstreamVLC(uint bufferSize = 50, EncoderParameters *p = nullptr);

    ~EncBitstreamVLC();


    uint write4DBlock(int *block, int size, const vector<LRE_struct> &lre, EncBitstreamVLC *prefix = nullptr);

    int writeGolomb_si(int val);

    int writeGolomb_ui(int val);

    Byte temp();

    int writeSyntaxElement(int info, int len);

    virtual EncBitstreamVLC &operator+=(const EncBitstreamVLC &temp);

    uint getSize() override;

    void reset() override;

    void printBitstream();

private:
    void writeUVLC2buffer(syntaxelement_enc *se);

    static int symbol2uvlc(syntaxelement_enc *sym);

    static void ue_linfo(int ue, int dummy, int *len, int *info);

    int write_ue_v(int value);

    int write_se_v(int value);

    int writeCoeff4x4_CAVLC_normal___H264(int *v, int start, uint size);

    int writeSyntaxElement_TotalZeros(syntaxelement_enc *se);

    int writeSyntaxElement_Level_VLCN(syntaxelement_enc *se, int vlc);

    int symbol2vlc(syntaxelement_enc *sym);

    int writeSyntaxElement_Run(syntaxelement_enc *se);

    int writeSyntaxElement_NumCoeffTrailingOnes(syntaxelement_enc *se);

    int writeSyntaxElement_Level_VLC1(syntaxelement_enc *se);

    int writeSyntaxElement_VLC(syntaxelement_enc *se);

    void se_linfo(int se, int dummy, int *len, int *info);

    void appendBitstream(const EncBitstreamVLC &temp);

};

#endif //LF_CODEC_ENCBITSTREAMVLC_H
