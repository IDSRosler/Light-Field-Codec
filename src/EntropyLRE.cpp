//
// Created by igor on 09/04/2022.
//

#include "EntropyLRE.h"

EntropyLRE::EntropyLRE() {}

std::vector<LRE_Struct> EntropyLRE::encodeLREVector(std::vector<bool> &v) {
    std::vector<LRE_Struct> lre;
    LRE_Struct lreValue;
    while (!v.empty()){
        lreValue.level = v.back() ? 1 : 0;
        lreValue.run = countSymbolRun(v, lreValue.level);
        lre.push_back(lreValue);
    }
    return lre;
}

int EntropyLRE::countSymbolRun(std::vector<bool> &v, int symbol){
    int run = 0;
    while (!v.empty() && v.back() == symbol){
        run++;
        v.pop_back();
    }
    return run;
}

#include "EntropyLRE.h"
