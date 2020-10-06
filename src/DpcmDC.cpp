
#include <cmath>
#include <algorithm>
#include "DpcmDC.h"

static int i_mean(const std::vector<int> &vet) {
    if (vet.empty()) return 0;

    double sum = 0;
    for (int i : vet) {
        sum += (int) i;
    }

    return (int) std::round(sum / vet.size());
}

static int i_median(std::vector<int> vet) {
    size_t size = vet.size();

    if (size == 0) return 0;

    std::sort(vet.begin(), vet.end());

    double val = (size % 2 == 0) ?
                 ((vet[size / 2 - 1] + vet[size / 2]) / 2) :
                 (vet[size / 2]);

    return (int) std::round(val);
}


ValueDPCM_DC::ValueDPCM_DC(int dc, bool available) : dc(dc), available(available) {}

std::ostream &operator<<(std::ostream &os, const ValueDPCM_DC &dcpm) {
    os << "available: " << dcpm.available << "\tdc: " << dcpm.dc;
    return os;
}

DpcmDC::DpcmDC(uint resol_x) : resol_x(resol_x) {
    this->init_references();
}

DpcmDC::~DpcmDC() {
    l_references.erase(l_references.begin(), l_references.end());
    l_references.clear();
}

int DpcmDC::get_reference(uint x, uint y) {
    int reference{0};

    ValueDPCM_DC r2 = *(this->l_references.begin() + 1); // above
    ValueDPCM_DC r3 = *(this->l_references.end() - 1); // left
    ValueDPCM_DC r1 = *this->l_references.begin(); // left above

    if (x == 0) r1.available = r3.available = false;
    if (y == 0) r1.available = r2.available = false;

    int num_available = r1.available + r2.available + r3.available;

    if (num_available == 3) {
        reference = (short) i_median({r1.dc, r2.dc, r3.dc});
    } else if (num_available == 2) {
        int dc1, dc2;

        if (r1.available && r2.available) {
            dc1 = r1.dc;
            dc2 = r2.dc;
        } else if (r2.available && r3.available) {
            dc1 = r2.dc;
            dc2 = r3.dc;
        } else {
            dc1 = r1.dc;
            dc2 = r3.dc;
        }

        reference = (short) i_mean({dc1, dc2});

    } else if (num_available == 1) {
        reference = (r1.available) ? (r1.dc) : ((r2.available) ? (r2.dc) : (r3.dc));
    }

    return reference;
}

void DpcmDC::update(int curr, bool available) {
    this->l_references.erase(this->l_references.begin());
    this->l_references.emplace_back(curr, available);
}

const std::vector<ValueDPCM_DC> &DpcmDC::getLReferences() const {
    return l_references;
}

void DpcmDC::init_references() {
    for (int i = 0; i < this->resol_x + 1; ++i) this->l_references.emplace_back();
}


