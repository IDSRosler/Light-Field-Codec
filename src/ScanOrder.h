
#ifndef _LFCODEC_SCANORDER_H_
#define _LFCODEC_SCANORDER_H_

#include <cstdint>
#include <unordered_map>
#include <algorithm>
#include <string>
#include "Point4D.h"
#include "utils.h"
#include "Typedef.h"

class ScanOrder {
  std::unordered_map<std::string, std::vector<float>> sums_map;
  static std::unordered_map<std::string, std::vector<index_t>> orders_map;
  static bool is_loaded;
  bool is_estimator = false;

  std::string make_key(const int type, const Point4D& shape) const;
  void read_file();
  void write_file();
  auto get_order(const std::string& key);
public:
  ScanOrder();
  explicit ScanOrder(bool is_estimator);
  ~ScanOrder();
  void process_block(int type, const float *block, const Point4D &shape,
                     const Point4D &stride);

  std::vector<index_t> get_order(int type, Point4D& shape);
};

#endif
