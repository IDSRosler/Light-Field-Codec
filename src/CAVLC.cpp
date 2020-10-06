/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CAVLC.cpp
 * Author: mateus
 *
 * Created on August 4, 2018, 9:35 PM
 */

#include "CAVLC.h"

CAVLC::CAVLC(bool is15x15x15x15) {
    this->bitstream = "";
    this->lre = new LRE(is15x15x15x15);
}

CAVLC::CAVLC(const CAVLC &orig) {
}

CAVLC::~CAVLC() {
    delete lre;
}

int CAVLC::writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se) {
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

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif

    return (se->len);
}

string CAVLC::writeCoeff4x4_CAVLC_normal___H264(int *v, int start, uint size) {

    this->bitstream = "";

    //// "DEFINES"
    vector<LRE_struct> czi = this->lre->encodeCZI(v, start, size);
    czi.push_back({0, 0});

//    for (int it = 0; it < 16; it++) {
//            cout << "lev: " << czi[it].level << "\trun: " << czi[it].run << endl;
//        }

    //    Slice* currSlice = currMB->p_Slice;
    //    VideoParameters *p_Vid = currSlice->p_Vid;
    int no_bits = 0;
    SyntaxElement se;
    //    DataPartition *dataPart;
    //    const int *partMap = assignSE2partition[currSlice->partition_mode];

    int k, level = 1, run = 0, vlcnum;
    int numcoeff = 0, lastcoeff = 0, numtrailingones = 0;
    int numones = 0, totzeros = 0, zerosleft, numcoef;
    int numcoeff_vlc = 0;
    int code, level_two_or_higher;
    int dptype = 0;
    //    int nnz, cdc = 0, cac = 0;
    int max_coeff_num = 0;
    //    int subblock_x, subblock_y;
    //    int *mb_bits_coeff = &currMB->bits.mb_y_coeff;
#if TRACE
    char type[15];
#endif

    static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768}; // maximum vlc = 6


    //    int* pLevel = NULL;
    //    int* pRun = NULL;
    //
    //    switch (block_type) {
    //        case LUMA:
    //            max_coeff_num = 16;
    //
    //            pLevel = currSlice->cofAC[b8][b4][0];
    //            pRun = currSlice->cofAC[b8][b4][1];
    //#if TRACE
    //            sprintf(type, "%s", "Luma");
    //#endif
    //            dptype = (is_intra(currMB)) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
    //            break;
    //
    //        case CHROMA_AC:
    //            max_coeff_num = 15;
    //            mb_bits_coeff = &currMB->bits.mb_uv_coeff;
    //            cac = 1;
    //
    //            pLevel = currSlice->cofAC[b8][b4][0];
    //            pRun = currSlice->cofAC[b8][b4][1];
    //#if TRACE
    //            sprintf(type, "%s", "ChrAC");
    //#endif
    //            dptype = (is_intra(currMB)) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
    //            break;
    //
    //        case CHROMA_DC:
    //            max_coeff_num = p_Vid->num_cdc_coeff;
    //            mb_bits_coeff = &currMB->bits.mb_uv_coeff;
    //            cdc = 1;
    //
    //            pLevel = currSlice->cofDC[param + 1][0];
    //            pRun = currSlice->cofDC[param + 1][1];
    //#if TRACE
    //            sprintf(type, "%s", "ChrDC");
    //#endif
    //            dptype = (is_intra(currMB)) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
    //            break;
    //
    //        case LUMA_INTRA16x16AC:
    //            max_coeff_num = 15;
    //
    //            pLevel = currSlice->cofAC[b8][b4][0];
    //            pRun = currSlice->cofAC[b8][b4][1];
    //#if TRACE
    //            sprintf(type, "%s", "Lum16AC");
    //#endif
    //            dptype = SE_LUM_AC_INTRA;
    //            break;
    //
    //        case LUMA_INTRA16x16DC:
    max_coeff_num = 16;
    //
    //            pLevel = currSlice->cofDC[0][0];
    //            pRun = currSlice->cofDC[0][1];
    //#if TRACE
    //            sprintf(type, "%s", "Lum16DC");
    //#endif
    dptype = 7; //SE_LUM_DC_INTRA;
    //            break;
    //
    //
    //        default:
    //            error("writeCoeff4x4_CAVLC: Invalid block type", 600);
    //            break;
    //    }


    //    dataPart = &(currSlice->partArr[partMap[dptype]]);

    //    for (k = 0; (k <= ((cdc) ? p_Vid->num_cdc_coeff : 16)) && level != 0; k++) {

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

    //    if (!cdc) {
    //        if (!cac) {
    //            // luma
    //            subblock_x = ((b8 & 0x1) == 0) ? (((b4 & 0x1) == 0) ? 0 : 1) : (((b4 & 0x1) == 0) ? 2 : 3);
    //            // horiz. position for coeff_count context
    //            subblock_y = (b8 < 2) ? ((b4 < 2) ? 0 : 1) : ((b4 < 2) ? 2 : 3);
    //            // vert.  position for coeff_count context
    //            nnz = predict_nnz(currMB, LUMA, subblock_x, subblock_y);
    //        } else {
    //            // chroma AC
    //            subblock_x = param >> 4;
    //            subblock_y = param & 15;
    //            nnz = predict_nnz_chroma(currMB, subblock_x, subblock_y);
    //        }
    //        p_Vid->nz_coeff [currMB->mbAddrX ][subblock_x][subblock_y] = numcoeff;
    //
    //        numcoeff_vlc = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));
    //    } else {
    //        // chroma DC (has its own VLC)
    //        // numcoeff_vlc not relevant
    numcoeff_vlc = 0;
    //
    //        subblock_x = param;
    //        subblock_y = param;
    //    }

    if (!numcoeff)
        return this->bitstream;

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

            // update VLC table
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

    return this->bitstream;
}

int CAVLC::symbol2vlc(SyntaxElement *sym) {
    string t = decToBinary_str_n(sym->inf, sym->len);
    this->bitstream += t;

    //    int info_len = sym->len;
    //
    //    // Convert info into a bitpattern int
    //    sym->bitpattern = 0;
    //
    //    // vlc coding
    //    while (--info_len >= 0) {
    //        sym->bitpattern <<= 1;
    //        sym->bitpattern |= (0x01 & (sym->inf >> info_len));
    //    }
    return 0;
}

int CAVLC::writeSyntaxElement_VLC(SyntaxElement *se) {
    se->inf = se->value1;
    se->len = se->value2;
    symbol2vlc(se);

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif

    return (se->len);
}

int CAVLC::writeSyntaxElement_Level_VLC1(SyntaxElement *se) {
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

        /* Assert to make sure that the code fits in the VLC */
        /* make sure that we are in High Profile to represent level_prefix > 15 */
        if (numPrefix > 0 && false) {
            //error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
            se->len = 0x0000FFFF; // This can be some other big number
            return (se->len);
        }

        se->len = 28 + (numPrefix << 1);
    }

    symbol2vlc(se);

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif
    //
    return (se->len);
}

int CAVLC::writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc) {
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

        /* Assert to make sure that the code fits in the VLC */
        /* make sure that we are in High Profile to represent level_prefix > 15 */
        if (numPrefix > 0 && false) {
            //error( "level_prefix must be <= 15 except in High Profile\n",  1000 );
            se->len = 0x0000FFFF; // This can be some other big number
            return (se->len);
        }
        se->len = 28 + (numPrefix << 1);
    }

    symbol2vlc(se);

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif

    return (se->len);
}

int CAVLC::writeSyntaxElement_TotalZeros(SyntaxElement *se) {
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

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif

    return (se->len);
}

int CAVLC::writeSyntaxElement_Run(SyntaxElement *se) {
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

    //    writeUVLC2buffer(se, dp->bitstream);
    //
    //    if (se->type != SE_HEADER)
    //        dp->bitstream->write_flag = 1;
    //
    //#if TRACE
    //    if (dp->bitstream->trace_enabled)
    //        trace2out(se);
    //#endif

    return (se->len);
}

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    buffer containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param bitcount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */
int CAVLC::ShowBits(Byte buffer[], int totbitoffset, int bitcount, int numbits) {
    if ((totbitoffset + numbits) > bitcount) {
        return -1;
    } else {
        int bitoffset = 7 - (totbitoffset & 0x07); // bit from start of byte
        int byteoffset = (totbitoffset >> 3); // byte from start of buffer
        Byte *curbyte = &(buffer[byteoffset]);
        int inf = 0;

        while (numbits--) {
            inf <<= 1;
            inf |= ((*curbyte) >> (bitoffset--)) & 0x01;

            if (bitoffset == -1) { //Move onto next byte to get all of numbits
                curbyte++;
                bitoffset = 7;
            }
        }
        return inf; // return absolute offset in bit from start of frame
    }
}

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer (Threshold based)
 *
 * \param inf
 *    bytes to extract numbits from with bitoffset already applied
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */
int CAVLC::ShowBitsThres(int inf, int numbits) {
    return ((inf) >> ((sizeof(Byte) * 24) - (numbits)));
    /*
    if ((numbits + 7) > bitcount)
    {
      return -1;
    }
    else
    {
      //Worst case scenario is that we will need to traverse 3 bytes
      inf >>= (sizeof(byte)*8)*3 - numbits;
    }

    return inf; //Will be a small unsigned integer so will not need any conversion when returning as int
     */
}

/*!
 ************************************************************************
 * \brief
 *    code from bitstream (2d tables)
 ************************************************************************
 */
int CAVLC::code_from_bitstream_2d(SyntaxElement *sym,
                                  Bitstream *currStream,
                                  const Byte *lentab,
                                  const Byte *codtab,
                                  int tabwidth,
                                  int tabheight,
                                  int *code) {
    int i, j;
    const Byte *len = &lentab[0], *cod = &codtab[0];

    int *frame_bitoffset = &currStream->frame_bitoffset;
    Byte *buf = &currStream->streamBuffer[*frame_bitoffset >> 3];

    //Apply bitoffset to three bytes (maximum that may be traversed by ShowBitsThres)
    unsigned int inf = ((*buf) << 16) + (*(buf + 1) << 8) + *(buf +
                                                              2); //Even at the end of a stream we will still be pulling out of allocated memory as alloc is done by MAX_CODED_FRAME_SIZE
    inf <<= (*frame_bitoffset & 0x07); //Offset is constant so apply before extracting different numbers of bits
    inf &= 0xFFFFFF; //Arithmetic shift so wipe any sign which may be extended inside ShowBitsThres

    // this VLC decoding method is not optimized for speed
    for (j = 0; j < tabheight; j++) {
        for (i = 0; i < tabwidth; i++) {
            if ((*len == 0) || (this->ShowBitsThres(inf, (int) *len) != *cod)) {
                ++len;
                ++cod;
            } else {
                sym->len = *len;
                *frame_bitoffset += *len; // move bitstream pointer
                *code = *cod;
                sym->value1 = i;
                sym->value2 = j;
                return 0; // found code and return
            }
        }
    }
    return -1; // failed to find code
}

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param info
 *    returns value of the read bits
 * \param bitcount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */
int CAVLC::GetBits(Byte buffer[], int totbitoffset, int *info, int bitcount,
                   int numbits) {
    if ((totbitoffset + numbits) > bitcount) {
        return -1;
    } else {
        int bitoffset = 7 - (totbitoffset & 0x07); // bit from start of byte
        int byteoffset = (totbitoffset >> 3); // byte from start of buffer
        int bitcounter = numbits;
        Byte *curbyte = &(buffer[byteoffset]);
        int inf = 0;

        while (numbits--) {
            inf <<= 1;
            inf |= ((*curbyte) >> (bitoffset--)) & 0x01;
            if (bitoffset == -1) { //Move onto next byte to get all of numbits
                curbyte++;
                bitoffset = 7;
            }
            // Above conditional could also be avoided using the following:
            // curbyte   -= (bitoffset >> 3);
            // bitoffset &= 0x07;
        }
        *info = inf;

        return bitcounter; // return absolute offset in bit from start of frame
    }
}

/*!
 ************************************************************************
 * \brief
 *    read  Run codeword from UVLC-partition
 ************************************************************************
 */
int CAVLC::readSyntaxElement_Run(SyntaxElement *sym, Bitstream *currStream) {
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
    int code;
    int vlcnum = sym->value1;
    int retval = code_from_bitstream_2d(sym, currStream, &lentab[vlcnum][0], &codtab[vlcnum][0], 16, 1, &code);

    if (retval) {
        printf("ERROR: failed to find Run\n");
        exit(-1);
    }

#if TRACE
    tracebits2(sym->tracestring, sym->len, code);
#endif

    return retval;
}

/*!
 ************************************************************************
 * \brief
 *    read Total Zeros codeword from UVLC-partition
 ************************************************************************
 */
int CAVLC::readSyntaxElement_TotalZeros(SyntaxElement *sym, Bitstream *currStream) {
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

    int code;
    int vlcnum = sym->value1;
    int retval = code_from_bitstream_2d(sym, currStream, &lentab[vlcnum][0], &codtab[vlcnum][0], 16, 1, &code);

    if (retval) {
        printf("ERROR: failed to find Total Zeros !cdc\n");
        exit(-1);
    }

#if TRACE
    tracebits2(sym->tracestring, sym->len, code);
#endif

    return retval;
}

/*!
 ************************************************************************
 * \brief
 *    read NumCoeff/TrailingOnes codeword from UVLC-partition
 ************************************************************************
 */

int CAVLC::readSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *sym,
                                                  Bitstream *currStream) {
    int frame_bitoffset = currStream->frame_bitoffset;
    int BitstreamLengthInBytes = currStream->bitstream_length;
    int BitstreamLengthInBits = (BitstreamLengthInBytes << 3) + 7;
    Byte *buf = currStream->streamBuffer;

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

    int retval = 0, code;
    int vlcnum = sym->value1;
    // vlcnum is the index of Table used to code coeff_token
    // vlcnum==3 means (8<=nC) which uses 6bit FLC

    if (vlcnum == 3) {
        // read 6 bit FLC
        //code = ShowBits(buf, frame_bitoffset, BitstreamLengthInBytes, 6);
        code = ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, 6);
        currStream->frame_bitoffset += 6;
        sym->value2 = (code & 3);
        sym->value1 = (code >> 2);

        if (!sym->value1 && sym->value2 == 3) {
            // #c = 0, #t1 = 3 =>  #c = 0
            sym->value2 = 0;
        } else
            sym->value1++;

        sym->len = 6;
    } else {
        //retval = code_from_bitstream_2d(sym, currStream, &lentab[vlcnum][0][0], &codtab[vlcnum][0][0], 17, 4, &code);
        retval = code_from_bitstream_2d(sym, currStream, lentab[vlcnum][0], codtab[vlcnum][0], 17, 4, &code);
        if (retval) {
            printf("ERROR: failed to find NumCoeff/TrailingOnes\n");
            assert(false);
        }
    }

#if TRACE
    snprintf(sym->tracestring, TRACESTRING_SIZE, "%s # c & tr.1s vlc=%d #c=%d #t1=%d",
            type, vlcnum, sym->value1, sym->value2);
    tracebits2(sym->tracestring, sym->len, code);
#endif

    return retval;
}

/*!
 ************************************************************************
 * \brief
 *    read FLC codeword from UVLC-partition
 ************************************************************************
 */
int CAVLC::readSyntaxElement_FLC(SyntaxElement *sym, Bitstream *currStream) {
    int BitstreamLengthInBits = (currStream->bitstream_length << 3) + 7;

    if ((GetBits(currStream->streamBuffer, currStream->frame_bitoffset, &(sym->inf), BitstreamLengthInBits, sym->len)) <
        0)
        return -1;

    sym->value1 = sym->inf;
    currStream->frame_bitoffset += sym->len; // move bitstream pointer

#if TRACE
    tracebits2(sym->tracestring, sym->len, sym->inf);
#endif

    return 1;
}

/*!
 ************************************************************************
 * \brief
 *    read Level VLC0 codeword from UVLC-partition
 ************************************************************************
 */
int CAVLC::readSyntaxElement_Level_VLC0(SyntaxElement *sym, Bitstream *currStream) {
    int frame_bitoffset = currStream->frame_bitoffset;
    int BitstreamLengthInBytes = currStream->bitstream_length;
    int BitstreamLengthInBits = (BitstreamLengthInBytes << 3) + 7;
    Byte *buf = currStream->streamBuffer;
    int len = 1, sign = 0, level = 0, code = 1;

    while (!ShowBits(buf, frame_bitoffset++, BitstreamLengthInBits, 1))
        len++;

    if (len < 15) {
        sign = (len - 1) & 1;
        level = ((len - 1) >> 1) + 1;
    } else if (len == 15) {
        // escape code
        code <<= 4;
        code |= ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, 4);
        len += 4;
        frame_bitoffset += 4;
        sign = (code & 0x01);
        level = ((code >> 1) & 0x07) + 8;
    } else if (len >= 16) {
        // escape code
        int addbit = (len - 16);
        int offset = (2048 << addbit) - 2032;
        len -= 4;
        code = ShowBits(buf, frame_bitoffset, BitstreamLengthInBits, len);
        sign = (code & 0x01);
        frame_bitoffset += len;
        level = (code >> 1) + offset;

        code |= (1 << (len)); // for display purpose only
        len += addbit + 16;
    }

    sym->inf = (sign) ? -level : level;
    sym->len = len;

#if TRACE
    tracebits2(sym->tracestring, sym->len, code);
#endif

    currStream->frame_bitoffset = frame_bitoffset;
    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    read Level VLC codeword from UVLC-partition
 ************************************************************************
 */
int CAVLC::readSyntaxElement_Level_VLCN(SyntaxElement *sym, int vlc, Bitstream *currStream) {
    int frame_bitoffset = currStream->frame_bitoffset;
    int BitstreamLengthInBytes = currStream->bitstream_length;
    int BitstreamLengthInBits = (BitstreamLengthInBytes << 3) + 7;
    Byte *buf = currStream->streamBuffer;

    int levabs, sign;
    int len = 1;
    int code = 1, sb;

    int shift = vlc - 1;

    // read pre zeros
    while (!ShowBits(buf, frame_bitoffset++, BitstreamLengthInBits, 1))
        len++;

    frame_bitoffset -= len;

    if (len < 16) {
        levabs = ((len - 1) << shift) + 1;

        // read (vlc-1) bits -> suffix
        if (shift) {
            sb = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, shift);
            code = (code << (shift)) | sb;
            levabs += sb;
            len += (shift);
        }

        // read 1 bit -> sign
        sign = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, 1);
        code = (code << 1) | sign;
        len++;
    } else // escape
    {
        int addbit = len - 5;
        int offset = (1 << addbit) + (15 << shift) - 2047;

        sb = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, addbit);
        code = (code << addbit) | sb;
        len += addbit;

        levabs = sb + offset;

        // read 1 bit -> sign
        sign = ShowBits(buf, frame_bitoffset + len, BitstreamLengthInBits, 1);

        code = (code << 1) | sign;

        len++;
    }

    sym->inf = (sign) ? -levabs : levabs;
    sym->len = len;

    currStream->frame_bitoffset = frame_bitoffset + len;

#if TRACE
    tracebits2(sym->tracestring, sym->len, code);
#endif

    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    Reads coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */
void CAVLC::read_coeff_4x4_CAVLC(Bitstream *currStream, int i, int j, int levarr[16], int runarr[16],
                                 int *number_coefficients) {
    //    Slice *currSlice = currMB->p_Slice;
    //    VideoParameters *p_Vid = currMB->p_Vid;
    //    int mb_nr = currMB->mbAddrX;
    SyntaxElement currSE;
    //    DataPartition *dP;
    //    const byte *partMap = assignSE2partition[currSlice->dp_mode];

    int k, code, vlcnum;
    int numcoeff = 0, numtrailingones;
    int level_two_or_higher;
    int numones, totzeros, abslevel, cdc = 0, cac = 0;
    int zerosleft, ntr, dptype = 0;
    int max_coeff_num = 0, nnz;
    char type[15];
    static const int incVlc[] = {0, 3, 6, 12, 24, 48, 32768}; // maximum vlc = 6

    /*switch (block_type) {
        case LUMA:
            max_coeff_num = 16;
            TRACE_PRINTF("Luma");
            dptype = (currMB->is_intra_block == TRUE) ? SE_LUM_AC_INTRA : SE_LUM_AC_INTER;
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
        case LUMA_INTRA16x16DC:
            max_coeff_num = 16;
            TRACE_PRINTF("Lum16DC");
            dptype = SE_LUM_DC_INTRA;
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
        case LUMA_INTRA16x16AC:
            max_coeff_num = 15;
            TRACE_PRINTF("Lum16AC");
            dptype = SE_LUM_AC_INTRA;
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
        case CHROMA_DC:
            max_coeff_num = p_Vid->num_cdc_coeff;
            cdc = 1;
            TRACE_PRINTF("ChrDC");
            dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
        case CHROMA_AC:
            max_coeff_num = 15;
            cac = 1;
            TRACE_PRINTF("ChrAC");
            dptype = (currMB->is_intra_block == TRUE) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
        default:
            error("read_coeff_4x4_CAVLC: invalid block type", 600);
            p_Vid->nz_coeff[mb_nr][0][j][i] = 0;
            break;
    }*/
    max_coeff_num = 16;
    dptype = 7; //SE_LUM_DC_INTRA;

    currSE.type = dptype;
    //    dP = &(currSlice->partArr[partMap[dptype]]);
    //    currStream = dP->bitstream;

    if (!cdc) {
        // luma or chroma AC
        nnz = 0; //(!cac) ? predict_nnz(currMB, LUMA, i << 2, j << 2) : predict_nnz_chroma(currMB, i, ((j - 4) << 2));

        currSE.value1 = (nnz < 2) ? 0 : ((nnz < 4) ? 1 : ((nnz < 8) ? 2 : 3));

        readSyntaxElement_NumCoeffTrailingOnes(&currSE, currStream);

        numcoeff = currSE.value1;
        numtrailingones = currSE.value2;

        //        p_Vid->nz_coeff[mb_nr][0][j][i] = (byte) numcoeff;
    }

    //    memset(levarr, 0, max_coeff_num * sizeof (int));
    //    memset(runarr, 0, max_coeff_num * sizeof (int));

    numones = numtrailingones;
    *number_coefficients = numcoeff;

    if (numcoeff) {
        if (numtrailingones) {
            currSE.len = numtrailingones;

#if TRACE
            snprintf(currSE.tracestring,
                    TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)", type, i, j);
#endif

            readSyntaxElement_FLC(&currSE, currStream);

            code = currSE.inf;
            ntr = numtrailingones;
            for (k = numcoeff - 1; k > numcoeff - 1 - numtrailingones; k--) {
                ntr--;
                levarr[k] = (code >> ntr) & 1 ? -1 : 1;
            }
        }

        // decode levels
        level_two_or_higher = (numcoeff > 3 && numtrailingones == 3) ? 0 : 1;
        vlcnum = (numcoeff > 10 && numtrailingones < 3) ? 1 : 0;

        for (k = numcoeff - 1 - numtrailingones; k >= 0; k--) {

#if TRACE
            snprintf(currSE.tracestring,
                    TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d ", type, i, j, k, vlcnum);
#endif

            if (vlcnum == 0)
                readSyntaxElement_Level_VLC0(&currSE, currStream);
            else
                readSyntaxElement_Level_VLCN(&currSE, vlcnum, currStream);

            if (level_two_or_higher) {
                currSE.inf += (currSE.inf > 0) ? 1 : -1;
                level_two_or_higher = 0;
            }

            levarr[k] = currSE.inf;
            abslevel = iabs(levarr[k]);
            if (abslevel == 1)
                ++numones;

            // update VLC table
            if (abslevel > incVlc[vlcnum])
                ++vlcnum;

            if (k == numcoeff - 1 - numtrailingones && abslevel > 3)
                vlcnum = 2;
        }

        if (numcoeff < max_coeff_num) {
            // decode total run
            vlcnum = numcoeff - 1;
            currSE.value1 = vlcnum;

#if TRACE
            snprintf(currSE.tracestring,
                    TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d ", type, i, j, vlcnum);
#endif
            if (false /*cdc*/)
                //                readSyntaxElement_TotalZerosChromaDC(p_Vid, &currSE, currStream);;
                int a = 0;
            else
                readSyntaxElement_TotalZeros(&currSE, currStream);

            totzeros = currSE.value1;
        } else {
            totzeros = 0;
        }

        // decode run before each coefficient
        zerosleft = totzeros;
        i = numcoeff - 1;

        if (zerosleft > 0 && i > 0) {
            do {
                // select VLC for runbefore
                vlcnum = imin(zerosleft - 1, RUNBEFORE_NUM_M1);

                currSE.value1 = vlcnum;
#if TRACE
                snprintf(currSE.tracestring,
                        TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d ",
                        type, i, j, i, vlcnum);
#endif

                readSyntaxElement_Run(&currSE, currStream);
                runarr[i] = currSE.value1;

                zerosleft -= runarr[i];
                i--;
            } while (zerosleft != 0 && i != 0);
        }
        runarr[i] = zerosleft;
    } // if numcoeff
}

/*!
 ************************************************************************
 * \brief
 *    mapping rule for ue(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    number in the code table
 ************************************************************************
 */
void CAVLC::linfo_ue(int len, int info, int *value1) {
    //assert ((len >> 1) < 32);
    *value1 = (int) (((unsigned int) 1 << (len >> 1)) + (unsigned int) (info) - 1);
}

/*!
 ************************************************************************
 * \brief
 *    mapping rule for se(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    signed mvd
 ************************************************************************
 */
void CAVLC::linfo_se(int len, int info, int *value1) {
    //assert ((len >> 1) < 32);
    unsigned int n = ((unsigned int) 1 << (len >> 1)) + (unsigned int) info - 1;
    *value1 = (n + 1) >> 1;
    if ((n & 0x01) == 0) // lsb is signed bit
        *value1 = -*value1;
}

/*!
 ************************************************************************
 * \brief
 *  read one exp-golomb VLC symbol
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param  info
 *    returns the value of the symbol
 * \param bytecount
 *    buffer length
 * \return
 *    bits read
 ************************************************************************
 */
int CAVLC::GetVLCSymbol(Byte buffer[], int totbitoffset, int *info, int bytecount) {
    long byteoffset = (totbitoffset >> 3); // byte from start of buffer
    int bitoffset = (7 - (totbitoffset & 0x07)); // bit from start of byte
    int bitcounter = 1;
    int len = 0;
    Byte *cur_byte = &(buffer[byteoffset]);
    int ctr_bit = ((*cur_byte) >> (bitoffset)) & 0x01; // control bit for current bit posision

    while (ctr_bit == 0) { // find leading 1 bit
        len++;
        bitcounter++;
        bitoffset--;
        bitoffset &= 0x07;
        cur_byte += (bitoffset == 7);
        byteoffset += (bitoffset == 7);
        ctr_bit = ((*cur_byte) >> (bitoffset)) & 0x01;
    }

    if (byteoffset + ((len + 7) >> 3) > bytecount)
        return -1;
    else {
        // make infoword
        int inf = 0; // shortest possible code is 1, then info is always 0

        while (len--) {
            bitoffset--;
            bitoffset &= 0x07;
            cur_byte += (bitoffset == 7);
            bitcounter++;
            inf <<= 1;
            inf |= ((*cur_byte) >> (bitoffset)) & 0x01;
        }

        *info = inf;
        return bitcounter; // return absolute offset in bit from start of frame
    }
}
