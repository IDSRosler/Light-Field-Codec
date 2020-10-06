#include "EncBitstreamVLC.h"

EncBitstreamVLC::EncBitstreamVLC(uint bufferSize, EncoderParameters *p) :
        EncBitstreamBuffer(bufferSize), size(p->dim_block.x) {

    this->lre = new LRE(p->dim_block.getNSamples() == 15 * 15 * 15 * 15);
}

int EncBitstreamVLC::writeGolomb_ui(int val) {
    return write_ue_v(val);
}

int EncBitstreamVLC::writeGolomb_si(int val) {
    return write_se_v(val);
}

EncBitstreamVLC::~EncBitstreamVLC() {
}

void EncBitstreamVLC::writeUVLC2buffer(syntaxelement_enc *se) {

    unsigned int mask = 1 << (se->len - 1);
    Byte *byte_buf = &this->byte_buf;
    uint *bits_to_go = &this->bits_to_go;
    int i;

    // Add the new bits to the bitstream.
    // Write out a byte if it is full
    if (se->len < 33) {
        for (i = 0; i < se->len; i++) {
            *byte_buf <<= 1u;

            if (se->bitpattern & mask)
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
        for (i = 0; i < (se->len - 32); i++) {
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

            if (se->bitpattern & mask)
                *byte_buf |= 1u;

            mask >>= 1u;

            if ((--(*bits_to_go)) == 0) {
                *bits_to_go = 8;
                this->buffer[this->byte_pos++] = *byte_buf;
                *byte_buf = 0;
            }
        }
    }
}

int EncBitstreamVLC::symbol2uvlc(syntaxelement_enc *sym) {
    uint suffix_len = sym->len >> 1u;
    //assert (suffix_len < 32);
    suffix_len = (1u << suffix_len);
    sym->bitpattern = suffix_len | (sym->inf & (suffix_len - 1));
    return 0;
}

void EncBitstreamVLC::ue_linfo(int ue, int dummy, int *len, int *info) {
    int i, nn = (ue + 1u) >> 1;

    for (i = 0; i < 33 && nn != 0; i++) {
        nn >>= 1u;
    }
    *len = (i << 1u) + 1;
    *info = ue + 1 - (1u << i);
}

void EncBitstreamVLC::se_linfo(int se, int dummy, int *len, int *info) {
    int sign = (se <= 0) ? 1 : 0;
    int n = iabs(se) << 1; //  n+1 is the number in the code table.  Based on this we find length and info
    int nn = (n >> 1);
    int i;
    for (i = 0; i < 33 && nn != 0; i++) {
        nn >>= 1;
    }
    *len = (i << 1) + 1;
    *info = n - (1 << i) + sign;
}

int EncBitstreamVLC::write_ue_v(int value) {
    syntaxelement_enc sym{};

    sym.value1 = value;
    sym.value2 = 0;

    ue_linfo(sym.value1, sym.value2, &(sym.len), &(sym.inf));
    symbol2uvlc(&sym);

    writeUVLC2buffer(&sym);

    return sym.len;
}

int EncBitstreamVLC::write_se_v(int value) {
    syntaxelement_enc sym{};
    sym.value1 = value;
    sym.value2 = 0;

    se_linfo(sym.value1, sym.value2, &(sym.len), &(sym.inf));
    symbol2uvlc(&sym);

    writeUVLC2buffer(&sym);

    return (sym.len);
}


int EncBitstreamVLC::writeCoeff4x4_CAVLC_normal___H264(int *v, int start, uint size) {

    //// "DEFINES"
    std::vector<LRE_struct> czi = this->lre->encodeCZI(v, start, size);
    czi.push_back({0, 0});

//    for (int it = 0; it < 16; it++) {
//            cout << "lev: " << czi[it].level << "\trun: " << czi[it].run << endl;
//        }

    int no_bits = 0;
    syntaxelement_enc se;

    int k, level = 1, run = 0, vlcnum;
    int numcoeff = 0, lastcoeff = 0, numtrailingones = 0;
    int numones = 0, totzeros = 0, zerosleft, numcoef;
    int numcoeff_vlc = 0;
    int code, level_two_or_higher;
    int dptype = 0;
    int max_coeff_num = 0;

#if TRACE
    char type[15];
#endif

    static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768}; // maximum vlc = 6

    max_coeff_num = 16;

    dptype = 7; //SE_LUM_DC_INTRA;


    for (k = 0; (k <= 16 && level != 0); k++) {
        level = czi[k].level; // level
        run = czi[k].run; // run

        if (level) {
            totzeros += run;
            if (iabs(level) == 1) {
                numones++;
                numtrailingones++;
                numtrailingones = imin(numtrailingones, 3); // clip to 3
            } else {
                numtrailingones = 0;
            }
            numcoeff++;
            lastcoeff = k;
        }
    }

    numcoeff_vlc = 0;

    if (!numcoeff)
        return 0;

    se.type = dptype;

    se.value1 = numcoeff;
    se.value2 = numtrailingones;
    se.len = numcoeff_vlc; /* use len to pass vlcnum */


    //    cout << "writeSyntaxElement_NumCoeffTrailingOnes: \n";
    writeSyntaxElement_NumCoeffTrailingOnes(&se);

    //    *mb_bits_coeff += se.len;
    no_bits += se.len;


    if (numcoeff) {
        code = 0;
        for (k = lastcoeff; k > lastcoeff - numtrailingones; k--) {
            //            level = pLevel[k]; // level
            level = czi[k].level; // level
            code <<= 1;

            code |= (level < 0);
        }

        if (numtrailingones) {
            se.type = dptype;

            se.value2 = numtrailingones;
            se.value1 = code;

#if TRACE
            snprintf(se.tracestring,
                    TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)",
                    type, subblock_x, subblock_y);
#endif
            //            cout << "writeSyntaxElement_VLC: \n";
            writeSyntaxElement_VLC(&se);
            //            *mb_bits_coeff += se.len;
            no_bits += se.len;

        }

        // encode levels
        level_two_or_higher = (numcoeff > 3 && numtrailingones == 3) ? 0 : 1;

        vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

        for (k = lastcoeff - numtrailingones; k >= 0; k--) {
            //            level = pLevel[k]; // level
            level = czi[k].level; // level

            se.value1 = level;
            se.type = dptype;

#if TRACE
            snprintf(se.tracestring,
                    TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d lev=%3d",
                    type, subblock_x, subblock_y, k, vlcnum, level);
#endif

            if (level_two_or_higher) {
                level_two_or_higher = 0;

                if (se.value1 > 0)
                    se.value1--;
                else
                    se.value1++;
            }

            //    encode level

            if (vlcnum == 0) {
                //                cout << "writeSyntaxElement_Level_VLC1: \n";

                writeSyntaxElement_Level_VLC1(&se);
            } else {
                //                cout << "writeSyntaxElement_Level_VLCN: \n";

                writeSyntaxElement_Level_VLCN(&se, vlcnum);
            }

            // update CAVLC table
            if (iabs(level) > incVlc[vlcnum])
                vlcnum++;

            if ((k == lastcoeff - numtrailingones) && iabs(level) > 3)
                vlcnum = 2;

            //            *mb_bits_coeff += se.len;
            no_bits += se.len;
        }

        // encode total zeroes
        if (numcoeff < max_coeff_num) {

            se.type = dptype;
            se.value1 = totzeros;

            vlcnum = numcoeff - 1;

            se.len = vlcnum;

#if TRACE
            snprintf(se.tracestring,
                    TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d totzeros=%3d",
                    type, subblock_x, subblock_y, vlcnum, totzeros);
#endif
            //            if (!cdc)
            //            cout << "writeSyntaxElement_TotalZeros: \n";
            writeSyntaxElement_TotalZeros(&se);
            //            else
            //                writeSyntaxElement_TotalZerosChromaDC(p_Vid, &se, dataPart);

            //            *mb_bits_coeff += se.len;
            no_bits += se.len;
        }

        // encode run before each coefficient
        zerosleft = totzeros;
        numcoef = numcoeff;
        for (k = lastcoeff; k >= 0; k--) {
            //            run = pRun[k]; //
            run = czi[k].run; // level

            se.value1 = run;
            se.type = dptype;

            // for last coeff, run is remaining totzeros
            // when zerosleft is zero, remaining coeffs have 0 run
            if ((!zerosleft) || (numcoeff <= 1))
                break;

            if (numcoef > 1 && zerosleft) {
                vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);
                se.len = vlcnum;

#if TRACE
                snprintf(se.tracestring,
                        TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d run=%2d",
                        type, subblock_x, subblock_y, k, vlcnum, run);
#endif
                //                cout << "writeSyntaxElement_Run: \n";
                writeSyntaxElement_Run(&se);

                //                *mb_bits_coeff += se.len;
                no_bits += se.len;

                zerosleft -= run;
                numcoef--;
            }
        }
    }
    return no_bits;
}

int EncBitstreamVLC::symbol2vlc(syntaxelement_enc *sym) {

    int info_len = sym->len;

    // Convert info into a bitpattern int
    sym->bitpattern = 0;

    // vlc coding
    while (--info_len >= 0) {
        sym->bitpattern <<= 1;
        sym->bitpattern |= (0x01 & (sym->inf >> info_len));
    }
    return 0;
}

int EncBitstreamVLC::writeSyntaxElement(int info, int len) {
    syntaxelement_enc se;

    se.value1 = info;
    se.value2 = len;

    se.inf = se.value1;
    se.len = se.value2;
    symbol2vlc(&se);

    writeUVLC2buffer(&se);

    return (se.len);
}

int EncBitstreamVLC::writeSyntaxElement_VLC(syntaxelement_enc *se) {
    se->inf = se->value1;
    se->len = se->value2;
    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

int EncBitstreamVLC::writeSyntaxElement_Run(syntaxelement_enc *se) {
    static const Byte lentab[TOTRUN_NUM][16] = {
            {1, 1},
            {1, 2, 2},
            {2, 2, 2, 2},
            {2, 2, 2, 3, 3},
            {2, 2, 3, 3, 3, 3},
            {2, 3, 3, 3, 3, 3, 3},
            {3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11},
    };

    static const Byte codtab[TOTRUN_NUM][16] = {
            {1, 0},
            {1, 1, 0},
            {3, 2, 1, 0},
            {3, 2, 1, 1, 0},
            {3, 2, 3, 2, 1, 0},
            {3, 0, 1, 3, 2, 5, 4},
            {7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    };
    int vlcnum = se->len;

    // se->value1 : run
    se->len = lentab[vlcnum][se->value1];
    se->inf = codtab[vlcnum][se->value1];

    if (se->len == 0) {
        printf("ERROR: (run) not valid: (%d)\n", se->value1);
        exit(-1);
    }

    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

int EncBitstreamVLC::writeSyntaxElement_NumCoeffTrailingOnes(syntaxelement_enc *se) {
    static const Byte lentab[3][4][17] = {
            { // 0702
                    {1, 6, 8, 9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16},
                    {0, 2, 6, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16},
                    {0, 0, 3, 7, 8, 9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16},
                    {0, 0, 0, 5, 6, 7, 8, 9, 10, 11, 13, 14, 14, 15, 15, 16, 16},
            },
            {
                    {2, 6, 6, 7, 8,  8,  9,  11, 11, 12, 12, 12, 13, 13, 13, 14, 14},
                    {0, 2, 5, 6, 6, 7,  8,  9,  11, 11, 12, 12, 13, 13, 14, 14, 14},
                    {0, 0, 3, 6, 6, 7, 8,  9,  11, 11, 12, 12, 13, 13, 13, 14, 14},
                    {0, 0, 0, 4, 4, 5, 6, 6, 7,  9,  11, 11, 12, 13, 13, 13, 14},
            },
            {
                    {4, 6, 6, 6, 7,  7,  7,  7,  8,  8,  9,  9,  9,  10, 10, 10, 10},
                    {0, 4, 5, 5, 5, 5,  6,  6,  7,  8,  8,  9,  9,  9,  10, 10, 10},
                    {0, 0, 4, 5, 5, 5, 6,  6,  7,  7,  8,  8,  9,  9,  10, 10, 10},
                    {0, 0, 0, 4, 4, 4, 4, 4, 5,  6,  7,  8,  8,  9,  10, 10, 10},
            },

    };

    static const Byte codtab[3][4][17] = {
            {
                    {1,  5,  7,  7, 7,  7,  15, 11, 8,  15, 11, 15, 11, 15, 11, 7, 4},
                    {0, 1,  4,  6,  6,  6, 6,  14, 10, 14, 10, 14, 10, 1,  14, 10, 6},
                    {0, 0, 1,  5,  5,  5, 5,  5, 13, 9,  13, 9, 13, 9, 13, 9,  5},
                    {0, 0, 0, 3,  3,  4,  4, 4, 4,  4,  12, 12, 8,  12, 8,  12, 8},
            },
            {
                    {3,  11, 7,  7, 7,  4,  7,  15, 11, 15, 11, 8,  15, 11, 7,  9, 7},
                    {0, 2,  7,  10, 6,  6, 6,  6,  14, 10, 14, 10, 14, 10, 11, 8,  6},
                    {0, 0, 3,  9,  5,  5, 5,  5, 13, 9,  13, 9, 13, 9, 6,  10, 5},
                    {0, 0, 0, 5,  4,  6,  8, 4, 4,  4,  12, 8,  12, 12, 8,  1,  4},
            },
            {
                    {15, 15, 11, 8, 15, 11, 9,  8,  15, 11, 15, 11, 8,  13, 9,  5, 1},
                    {0, 14, 15, 12, 10, 8, 14, 10, 14, 14, 10, 14, 10, 7,  12, 8,  4},
                    {0, 0, 13, 14, 11, 9, 13, 9, 13, 10, 13, 9, 13, 9, 11, 7,  3},
                    {0, 0, 0, 12, 11, 10, 9, 8, 13, 12, 12, 12, 8,  12, 10, 6,  2},
            },
    };
    int vlcnum = se->len;

    // se->value1 : numcoeff
    // se->value2 : numtrailingones

    if (vlcnum == 3) {
        se->len = 6; // 4 + 2 bit FLC
        if (se->value1 > 0) {
            se->inf = ((se->value1 - 1) << 2) | se->value2;
        } else {
            se->inf = 3;
        }
    } else {
        se->len = lentab[vlcnum][se->value2][se->value1];
        se->inf = codtab[vlcnum][se->value2][se->value1];
    }

    if (se->len == 0) {
        printf("ERROR: (numcoeff,trailingones) not valid: vlc=%d (%d, %d)\n",
               vlcnum, se->value1, se->value2);
        exit(-1);
    }

    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

int EncBitstreamVLC::writeSyntaxElement_Level_VLC1(syntaxelement_enc *se) {
    int level = se->value1;
    int sign = (level < 0 ? 1 : 0);
    int levabs = iabs(level);

    if (levabs < 8) {
        se->len = levabs * 2 + sign - 1;
        se->inf = 1;
    } else if (levabs < 16) {
        // escape code1
        se->len = 19;
        se->inf = 16 | ((levabs << 1) - 16) | sign;
    } else {
        int iMask = 4096, numPrefix = 0;
        int levabsm16 = levabs + 2032;

        // escape code2
        if ((levabsm16) >= 4096) {
            numPrefix++;
            while ((levabsm16) >= (4096 << numPrefix)) {
                numPrefix++;
            }
        }

        iMask <<= numPrefix;
        se->inf = iMask | ((levabsm16 << 1) - iMask) | sign;

        se->len = 28 + (numPrefix << 1);
    }

    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

int EncBitstreamVLC::writeSyntaxElement_Level_VLCN(syntaxelement_enc *se, int vlc) {
    int level = se->value1;
    int sign = (level < 0 ? 1 : 0);
    int levabs = iabs(level) - 1;

    int shift = vlc - 1;
    int escape = (15 << shift);

    if (levabs < escape) {
        int sufmask = ~((0xffffffff) << shift);
        int suffix = (levabs) & sufmask;

        se->len = ((levabs) >> shift) + 1 + vlc;
        se->inf = (2 << shift) | (suffix << 1) | sign;
    } else {
        int iMask = 4096;
        int levabsesc = levabs - escape + 2048;
        int numPrefix = 0;

        if ((levabsesc) >= 4096) {
            numPrefix++;
            while ((levabsesc) >= (4096 << numPrefix)) {
                numPrefix++;
            }
        }

        iMask <<= numPrefix;
        se->inf = iMask | ((levabsesc << 1) - iMask) | sign;

        se->len = 28 + (numPrefix << 1);
    }

    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

int EncBitstreamVLC::writeSyntaxElement_TotalZeros(syntaxelement_enc *se) {
    static const Byte lentab[TOTRUN_NUM][16] = {
            {1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9},
            {3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6},
            {4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6},
            {5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5},
            {4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5},
            {6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6},
            {6, 5, 3, 3, 3, 2, 3, 4, 3, 6},
            {6, 4, 5, 3, 2, 2, 3, 3, 6},
            {6, 6, 4, 2, 2, 3, 2, 5},
            {5, 5, 3, 2, 2, 2, 4},
            {4, 4, 3, 3, 1, 3},
            {4, 4, 2, 1, 3},
            {3, 3, 1, 2},
            {2, 2, 1},
            {1, 1},
    };

    static const Byte codtab[TOTRUN_NUM][16] = {
            {1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1},
            {7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0},
            {5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0},
            {3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0},
            {5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0},
            {1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0},
            {1, 1, 5, 4, 3, 3, 2, 1, 1, 0},
            {1, 1, 1, 3, 3, 2, 2, 1, 0},
            {1, 0, 1, 3, 2, 1, 1, 1,},
            {1, 0, 1, 3, 2, 1, 1,},
            {0, 1, 1, 2, 1, 3},
            {0, 1, 1, 1, 1},
            {0, 1, 1, 1},
            {0, 1, 1},
            {0, 1},
    };
    int vlcnum = se->len;

    // se->value1 : TotalZeros
    se->len = lentab[vlcnum][se->value1];
    se->inf = codtab[vlcnum][se->value1];

    if (se->len == 0) {
        assert(se->len == 0);
        printf("ERROR: (TotalZeros) not valid: (%d)\n", se->value1);
        exit(-1);
    }

    symbol2vlc(se);

    writeUVLC2buffer(se);

    return (se->len);
}

uint EncBitstreamVLC::write4DBlock(int *block, int size, const vector<LRE_struct> &lre, EncBitstreamVLC *prefix) {

    bool has_bitstream;
    uint bitstream_size = 0;

    if (this->entropy_type == 1) {
        //vector<LRE_struct> czi_t = this->lre->encodeCZI(block, 0, size);
        //bitstream_size += writeGolomb_ui(czi_t.size());

        for (auto &i : lre) {
            bitstream_size += writeGolomb_si(i.level) + writeGolomb_ui(i.run);
        }
        bitstream_size += writeGolomb_ui(0) + writeGolomb_ui(0);
    } else {
        assert(prefix != nullptr);

        for (int i = 0; i < size; i += 16) {
            has_bitstream = this->writeCoeff4x4_CAVLC_normal___H264(block, i, i + 16);
            prefix->writeSyntaxElement(has_bitstream, 1);
        }

        bitstream_size = prefix->getSize() + this->getSize();
    }

    return bitstream_size;
}


void EncBitstreamVLC::appendBitstream(const EncBitstreamVLC &temp) {
    for (int i = 0; i < temp.byte_pos; i++) {
        this->writeSyntaxElement(temp.buffer[i], 8);
    }
    this->writeSyntaxElement(temp.byte_buf, 8 - temp.bits_to_go);
}

void EncBitstreamVLC::printBitstream() {
    for (int i = 0; i < this->byte_pos; i++) {
        cout << bitset<8>(this->buffer[i]);
    }
    cout << endl;
}

Byte EncBitstreamVLC::temp() {
    return this->byte_buf;
}

EncBitstreamVLC &EncBitstreamVLC::operator+=(const EncBitstreamVLC &temp) {
    this->appendBitstream(temp);
    return *this;
}

uint EncBitstreamVLC::getSize() {
    return EncBitstreamBuffer::getSize();
}

void EncBitstreamVLC::reset() {
    EncBitstreamBuffer::reset();
}
