#ifndef LF_CODEC_LRE_H
#define LF_CODEC_LRE_H

#include "Point4D.h"
#include "Typedef.h"
#include <iostream>
#include <vector>
#include <unordered_map>

struct LRE_struct {
    int level;
    std::size_t run;
};

class LRE {
public:
    LRE() = default;
    LRE(bool is15x15x15x15);
    LRE(const Point4D& shape);
    LRE(const LRE &orig);

    virtual ~LRE();

    std::vector<LRE_struct> encodeLRE(const int *v, std::size_t end);
    std::vector<LRE_struct> encodeLRE(const int *v, std::size_t end,  std::vector<index_t> scan_order);
    void decodeLRE(const std::vector<LRE_struct>& v_lre, std::size_t size, int *into);
    void decodeLRE(const std::vector<LRE_struct>& v_lre, std::size_t size, int *into, const std::vector<index_t>& scan_order);

    std::vector<LRE_struct> encodeSequence(const int *v, std::size_t end);

    std::vector<LRE_struct> encodeCZI(const int *v, int start, std::size_t end);
    std::vector<LRE_struct> encodeCZI(const int *v, int start, std::size_t end, std::vector<index_t> scan_order);

    std::vector<index_t> getScanOrder() const;
private:
    std::vector<index_t> scan_order;

    int find_first_not_of(const int *v, std::size_t pos_orig, std::size_t end);
    int find_first_not_of(const int *v, std::size_t pos_orig, std::size_t end,  std::vector<index_t> scan_order);
    int countZeroBefore(const int *v, std::size_t pos, std::size_t end);
    int countZeroBefore(const int *v, std::size_t pos, std::size_t end, std::vector<index_t> scan_order);

};
#endif //LF_CODEC_LRE_H
