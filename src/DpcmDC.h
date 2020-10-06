
#ifndef DPCMDC_H
#define DPCMDC_H

#include <vector>
#include <ostream>

class ValueDPCM_DC {
public:
    bool available{false};
    int dc{0};

    ValueDPCM_DC() = default;

    ValueDPCM_DC(int dc, bool available);

    friend std::ostream &operator<<(std::ostream &os, const ValueDPCM_DC &dcpm);
};

class DpcmDC {
private:
    uint resol_x;
    std::vector<ValueDPCM_DC> l_references;

    void init_references();

public:
    DpcmDC(uint resol_x);

    ~DpcmDC();

    int get_reference(uint x, uint y);

    void update(int curr, bool available);

    const std::vector<ValueDPCM_DC> &getLReferences() const;

};

#endif //DPCMDC_H
