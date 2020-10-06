
#include "ScanOrder.h"
#include "Transform.h"

std::unordered_map<std::string, std::vector<index_t>> ScanOrder::orders_map;
bool ScanOrder::is_loaded;

ScanOrder::ScanOrder() : ScanOrder(false) {}

ScanOrder::ScanOrder(bool is_estimator) : is_estimator(is_estimator) {
  if (!is_estimator)
    read_file();
}


std::string ScanOrder::make_key(const int type, const Point4D& shape) const
{
  std::stringstream ss;
  auto t = (Transform::TransformType)type;
  switch (t) {
  case Transform::DST_I:
    ss << "DST_I(";
    break;
  case Transform::DST_VII:
    ss << "DST_VII(";
    break;
  case Transform::DCT_II:
  default:
    ss << "DCT_II(";
    break;
  }
  ss << shape.x << ";" << shape.y << ";" << shape.u << ";" << shape.v << ")";
  return ss.str();
}

void ScanOrder::process_block(int type, const float *block,
                              const Point4D &shape, const Point4D &stride) {
  if (!is_estimator)
      return;
  auto key = make_key(type, shape);
  bool is_block_new = false;
  std::vector<float> _block;
  auto size = stride.v * shape.v;



  try {
    _block = sums_map.at(key);
  } catch (...) {
    _block.resize(size);
  }

  FOREACH_4D_IDX(i, shape, stride) {
    auto old_value = is_block_new ? 0 : _block[i];
    _block[i] = (block[i] * block[i]) + old_value;
  }
  sums_map[key] = _block;
}

auto ScanOrder::get_order(const std::string& key) {
  std::vector acc_block = sums_map.at(key);
  std::vector<std::pair<float, index_t>> avg_block;
  std::vector<index_t> scan_order;
  auto size = acc_block.size();

  avg_block.resize(size);
  scan_order.resize(size);

  for (std::size_t i = 0; i < size; i++) {
    avg_block[i].first = acc_block[i];
    avg_block[i].second = i;
  }

  std::sort(avg_block.begin(), avg_block.end(),
            [](auto &a, auto &b) -> bool { return a.first > b.first; });

  for (std::size_t i = 0; i < size; i++)
    scan_order[i] = avg_block[i].second;

  orders_map[key] = scan_order;
  return scan_order;
}

std::vector<index_t> ScanOrder::get_order(int type, Point4D &shape) {
  auto key = make_key(type, shape);
  try {
    return orders_map.at(key);
  } catch (...) {
    return get_order(key);
  }
}

ScanOrder::~ScanOrder()
{
  if (is_estimator)
    write_file();
}

void ScanOrder::write_file() {
  std::vector<std::string> keys;
  std::vector<std::vector<index_t>> orders;

  std::transform(sums_map.begin(),
                 sums_map.end(),
                 std::back_inserter(keys),
                 [](auto &k){return k.first;});

  for (auto& it : keys) {
    auto order = get_order(it);
    orders.push_back(order);
  }

  std::string sep = ",";

  std::ofstream file;
  file.open("scan_order.csv");
  file << "Index";
  for (const auto& it : keys)
    file << sep << it;
  file << std::endl;

  std::size_t i = 0;
  std::size_t nop;
  do {
    std::stringstream ss;
    nop = 0;
    ss << i;
    for (auto & order : orders) {
      ss << sep;
      if (i < order.size())
        ss << order[i];
      else
        nop++;
    }
    ss << std::endl;
    if (nop < orders.size())
      file << ss.str();
    i++;
  } while (nop < orders.size());

  file.close();
}

void ScanOrder::read_file() {
  if (is_loaded)
    return;

  std::ifstream file("default_scan_order.csv");
  std::vector<std::pair<std::string, std::vector<index_t>>> result;

  if(!file.is_open()) throw std::runtime_error("Could not open file");

  std::string line, colname;

  if(file.good())
  {
    std::getline(file, line);

    std::stringstream ss(line);

    while(std::getline(ss, colname, ',')){
      result.push_back({colname, std::vector<index_t> {}});
    }
  }

  while(std::getline(file, line))
  {
    std::stringstream ss(line);
    int colIdx = 0;
    std::stringstream conv;
    int val;
    std::string str;

    while(std::getline(ss, str, ',')) {
      if (str != "") {
        conv << str;
        conv >> val;
        result.at(colIdx).second.push_back(val);
        conv.clear();
      }
      colIdx++;
    }
  }

  // Close file
  file.close();

  for (auto [colname, vec] : result) {
    orders_map[colname] = vec;
  }
  is_loaded = true;
}