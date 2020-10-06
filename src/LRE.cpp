#include "LRE.h"

LRE::LRE(bool is15x15x15x15) {
    this->scan_order = (is15x15x15x15) ? this->basic_scan_order_15x15x15x15 : this->basic_scan_order_15x15x13x13;
}

LRE::LRE(const LRE &orig) {
}

LRE::~LRE() {
}

std::vector<LRE_struct> LRE::encodeLRE(const int *v, uint end) {
    std::vector<LRE_struct> v_lre;
    uint n_ocorr = 0;

    for (uint i = 0; i < end; i += (n_ocorr + 1)) {
        n_ocorr = this->find_first_not_of(v, i, end);
        v_lre.push_back({v[this->scan_order[i]], n_ocorr + 1});
    }

    return v_lre;
}

std::vector<LRE_struct> LRE::encodeCZI(const int *v, int start, uint end) {
    std::vector<LRE_struct> v_lre;
    uint n_zeroBefore = 0;

    for (uint i = start; i < end; i += (n_zeroBefore + 1)) {
        n_zeroBefore = this->countZeroBefore(v, i, end);

        if (i + n_zeroBefore < end) {
            v_lre.push_back({v[this->scan_order[i + n_zeroBefore]], n_zeroBefore});
        }
    }

    return v_lre;
}

int LRE::find_first_not_of(const int *v, uint pos_orig, uint end) {
    uint ocorr = 0;
    uint pos = pos_orig + 1;

    while ((pos < end) and (v[this->scan_order[pos]] == v[this->scan_order[pos_orig]])) {
        ++ocorr;
        ++pos;
    }

    return ocorr;
}

int LRE::countZeroBefore(const int *v, uint pos, uint end) {
    uint zero = 0;
    while ((pos < end) and (v[this->scan_order[pos++]] == 0)) zero++;
    return zero;
}

const ushort *LRE::getScanOrder() const {
    return this->scan_order;
}
