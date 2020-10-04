#include "LRE.h"
#include "utils.h"
#include "EncoderParameters.h"

LRE::LRE(bool is15x15x15x15) {
  Point4D shape, stride;
  if (is15x15x15x15)
    shape.set(15, 15, 15, 15);
  else
    shape.set(15, 15, 13, 13);
  shape.updateNSamples();
  stride = make_stride(shape);
  // switch(EncoderParameters::parameters.scan_order)
  // {
  //   case EncoderParameters::Z_ORDER_CURVE:
  //     scan_order = generate_z_order_curve(shape, stride);
  //     break;
  //   case EncoderParameters::ZIG_ZAG_SCAN_ORDER:
  //     scan_order = generate_scan_order(shape, stride);
  //     break;
  //   case EncoderParameters::CUSTOM_SCAN_ORDER:
  //     FOREACH_4D_IDX(i, shape, stride)
  //       scan_order.push_back(i);
  // }
  scan_order = generate_scan_order(shape, stride);
}


LRE::LRE(const Point4D &shape) {
    Point4D stride = make_stride(shape);
    scan_order = generate_scan_order(shape, stride);
}

LRE::LRE(const LRE &orig) {
}

LRE::~LRE() {
}

std::vector<LRE_struct> LRE::encodeLRE(const int *v, std::size_t end) {
  return encodeLRE(v, end, scan_order);
}
std::vector<LRE_struct> LRE::encodeLRE(const int *v, std::size_t end, std::vector<index_t> scan_order) {
    std::vector<LRE_struct> v_lre;
    std::size_t n_ocorr = 0;

    for (std::size_t i = 0; i < end; i += (n_ocorr + 1)) {
        n_ocorr = this->find_first_not_of(v, i, end, scan_order);
        v_lre.push_back({v[scan_order[i]], n_ocorr + 1});
    }

    return v_lre;
}



void LRE::decodeLRE(const std::vector<LRE_struct>& v_lre, std::size_t size, int *into)
{
    decodeLRE(v_lre, size, into, scan_order);
}
void LRE::decodeLRE(const std::vector<LRE_struct>& v_lre, std::size_t size, int *into, const std::vector<index_t>& scan_order) {
    std::size_t run;
    std::size_t scan_order_pos = 0;
    std::size_t pos;

    for (const auto& obj: v_lre) {
        run = obj.run;
        while (run--) {
            do {
                pos = scan_order[scan_order_pos++];
            } while (pos >= size);
            into[pos] = obj.level;
        }
    }
}

std::vector<LRE_struct> LRE::encodeCZI(const int *v, int start, std::size_t end) {
    return encodeCZI(v, start, end, scan_order);
}


std::vector<LRE_struct> LRE::encodeCZI(const int *v, int start, std::size_t end, std::vector<index_t> scan_order)
{
    std::vector<LRE_struct> v_lre;
    std::size_t n_zeroBefore = 0;

    for (std::size_t i = start; i < end; i += (n_zeroBefore + 1)) {
        n_zeroBefore = countZeroBefore(v, i, end, scan_order);

        if (i + n_zeroBefore < end) {
            v_lre.push_back({v[scan_order[i + n_zeroBefore]], n_zeroBefore});
        }
    }

    return v_lre;
}

int LRE::find_first_not_of(const int *v, std::size_t pos_orig, std::size_t end) {
  return find_first_not_of(v, pos_orig, end, scan_order);
}

int LRE::find_first_not_of(const int *v, std::size_t pos_orig, std::size_t end, std::vector<index_t> scan_order) {
    std::size_t ocorr = 0;
    std::size_t pos = pos_orig + 1;

    while ((pos < end) && (v[scan_order[pos]] == v[scan_order[pos_orig]])) {
        ++ocorr;
        ++pos;
    }

    return ocorr;
}

int LRE::countZeroBefore(const int *v, std::size_t pos, std::size_t end) {
    return countZeroBefore(v, pos, end, scan_order);
}

int LRE::countZeroBefore(const int *v, std::size_t pos, std::size_t end, std::vector<index_t> scan_order) {
    std::size_t zero = 0;
    while ((pos < end) && (v[scan_order[pos++]] == 0)) zero++;
    return zero;
}

std::vector<index_t> LRE::getScanOrder() const {
    return this->scan_order;
}


