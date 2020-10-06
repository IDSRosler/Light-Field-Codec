#ifndef CAVLC_H
#define CAVLC_H

#include <limits.h>
#include <vector>
#include <iostream>
#include <assert.h>
#include "LRE.h"

typedef unsigned char Byte; //!< byte type definition
using namespace std;

#define TOTRUN_NUM       15
#define RUNBEFORE_NUM_M1  6

typedef struct syntaxelement_dec {
    int type; //!< type of syntax element for data part.
    int value1; //!< numerical value of syntax element
    int value2; //!< for blocked symbols, e.g. run/level
    int len; //!< length of code
    int inf; //!< info part of CAVLC code
    unsigned int bitpattern; //!< CAVLC bitpattern
    int context; //!< CABAC context
    int k; //!< CABAC context for coeff_count,uv

    //! for mapping of CAVLC to syntaxElement
    void (*mapping)(int len, int info, int *value1, int *value2);
} SyntaxElement;

static inline string decToBinary_str_n(int n, int sz) {
    string val = "";
    int binaryNum[1000];

    int i = 0;
    while (n > 0) {
        binaryNum[i] = n % 2;
        n = n / 2;
        ++i;
    }

    for (int j = i - 1; j >= 0; --j) {
        val += (binaryNum[j]) ? "1" : "0";
    }

    for (int i = val.size(); i < sz; ++i) {
        val = "0" + val;
    }

    return val;
}

static inline int iabs(int x) {
    static const int INT_BITS = (sizeof(int) * CHAR_BIT) - 1;
    int y = x >> INT_BITS;
    return (x ^ y) - y;
}

static inline int imin(int a, int b) {
    return ((a) < (b)) ? (a) : (b);
}

struct bit_stream_dec {
    // CABAC Decoding
    int read_len; //!< actual position in the codebuffer, CABAC only
    int code_len; //!< overall codebuffer length, CABAC only
    // CAVLC Decoding
    int frame_bitoffset; //!< actual position in the codebuffer, bit-oriented, CAVLC only
    int bitstream_length; //!< over codebuffer lnegth, byte oriented, CAVLC only
    // ErrorConcealment
    Byte *streamBuffer; //!< actual codebuffer for read bytes
    int ei_flag; //!< error indication, 0: no error, else unspecified error

    bit_stream_dec() {
        this->read_len = this->code_len = 0;
        this->frame_bitoffset = this->bitstream_length = 0;
        this->streamBuffer = nullptr;
        this->ei_flag = 0;
    };
};

typedef struct bit_stream_dec Bitstream;

class CAVLC {
public:
    string bitstream;

    CAVLC(bool is15x15x15x15);

    CAVLC(const CAVLC &orig);

    virtual ~CAVLC();

    int symbol2vlc(SyntaxElement *sym);

    string writeCoeff4x4_CAVLC_normal___H264(int *v, int start, uint size);

    int writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se);

    int writeSyntaxElement_VLC(SyntaxElement *se);

    int writeSyntaxElement_Level_VLC1(SyntaxElement *se);

    int writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc);

    int writeSyntaxElement_TotalZeros(SyntaxElement *se);

    int writeSyntaxElement_Run(SyntaxElement *se);

    // decoder
    int ShowBits(Byte buffer[], int totbitoffset, int bitcount, int numbits);

    int ShowBitsThres(int inf, int numbits);

    int code_from_bitstream_2d(SyntaxElement *sym,
                               Bitstream *currStream, const Byte *lentab, const Byte *codtab,
                               int tabwidth, int tabheight, int *code);

    int GetBits(Byte buffer[], int totbitoffset, int *info, int bitcount, int numbits);

    int readSyntaxElement_Run(SyntaxElement *sym, Bitstream *currStream);

    int readSyntaxElement_TotalZeros(SyntaxElement *sym, Bitstream *currStream);

    int readSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *sym, Bitstream *currStream);

    int readSyntaxElement_FLC(SyntaxElement *sym, Bitstream *currStream);

    int readSyntaxElement_Level_VLC0(SyntaxElement *sym, Bitstream *currStream);

    int readSyntaxElement_Level_VLCN(SyntaxElement *sym, int vlc, Bitstream *currStream);

    void
    read_coeff_4x4_CAVLC(Bitstream *currStream, int i, int j, int levarr[16], int runarr[16], int *number_coefficients);

    void linfo_ue(int len, int info, int *value1);

    void linfo_se(int len, int info, int *value1);

    int GetVLCSymbol(Byte buffer[], int totbitoffset, int *info, int bytecount);

private:
    LRE *lre;

};

#endif /* CAVLC_H */

