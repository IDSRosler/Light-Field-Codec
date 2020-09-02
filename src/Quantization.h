
#ifndef QUANTIZATION_H
#define QUANTIZATION_H

#include <cmath>
#include <tuple>
#include <unordered_map>
#include <functional>
#include "Typedef.h"
#include "Point4D.h"
#include "EncoderParameters.h"

class Quantization {

public:
  enum QuantizationVolumeType {
    HAIYAN = 1,
    LEE = 2,
    EXPONENTIAL = 3
  };

  Quantization() = default;

  Quantization(const Point4D &dimBlock, const EncoderParameters codec_parameters);

  virtual ~Quantization();

  void delete_volume();

  void set_lee_parameters(float C, float ai, float a0, float bi, float b0);

  void inverse(const float *input, float *output);

  void foward(const float *input, float *output);

  void inverse(const QuantizationVolumeType type, const float *input, float *output);

  void foward(const QuantizationVolumeType type, const float *input, float *output);

  float qp = 1;
  float *get_volume(const QuantizationVolumeType type);
private:
  float *volume_hayan = nullptr;
  float *volume_lee = nullptr;
  float lee_c = 10;
  float lee_ai = 1023;
  float lee_a0 = 1023;
  float lee_bi = 0.04;
  float lee_b0 = 0.01;
  Point4D dim_block;
  Point4D weight_100;
  
  inline static std::unordered_map<const Point4D, float *, Point4DHasher> cache_lee;
  inline static std::unordered_map<const Point4D, float *, Point4DHasher> cache_hayan;


  static float clip_min(float val, float min);

  float *generate_volume(const QuantizationVolumeType type);

};


#endif //QUANTIZATION_H
