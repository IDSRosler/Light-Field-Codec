//
// Created by igor on 09/04/2022.
//

#ifndef LF_CODEC_ENTROPYLRE_H
#define LF_CODEC_ENTROPYLRE_H

#include <vector>

#include "ArithmeticStructures.h"

class EntropyLRE {
public:
    EntropyLRE();
    std::vector<LRE_Struct> encodeLREVector(std::vector<bool> &v);

private:
    int countSymbolRun(std::vector<bool> &v, int symbol);
};


#endif //LF_CODEC_ENTROPYLRE_H
