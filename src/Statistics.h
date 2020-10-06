#ifndef STATISTICS_H
#define STATISTICS_H

#include <cmath>
#include <climits>
#include <vector>
#include <algorithm>
#include <tuple>
#include <fstream>
#include <cassert>
#include <cfloat>

#include "Point4D.h"
#include "LRE.h"


class Statistics {
public:

    explicit Statistics(const std::string &file = "");

    virtual ~Statistics();

    void write(Point4D &pos, Point4D &dimBlock, uint it_channel, std::vector<LRE_struct> &lre_result, uint bits_per_4D_Block);

    void write(Point4D &pos, Point4D &dimBlock, uint it_channel);

    double get_mean() const;

    double get_median() const;

    double get_variance() const;

    double get_std() const;

    unsigned long get_zeros() const;

    unsigned long get_ones() const;

    double get_entropy() const;

    double get_power() const;

    void compute(const std::vector<float> &input, const ushort *scan_order);

    void compute_sse(float *orig, float *ref, const Point4D &dim_block, const Point4D &stride_block);

    static float min_vet(const std::vector<float> &input);

    static float max_vet(const std::vector<float> &input);

    static float minAbs_vet(const std::vector<float> &input);

    static float maxAbs_vet(const std::vector<float> &input);

private:
    std::string sep{"\t"};

    std::ofstream file_out;
    unsigned long num_zeros{}, num_ones{};
    double v_mean{}, v_median{}, v_variance{}, v_std{};
    double v_entropy{}, v_energy{};
    float dc, sse{0}, v_max, v_min, v_maxAbs, v_minAbs;

    double run_mean, run_var, run_std;
    int posSO_last_nzero, posSO_last_nzeroone;
    int pos_last_nzero, pos_last_nzeroone; // so = scan order

    static double variance(std::vector<float> const &vet, double mean);

    static double median(std::vector<float> vet);

    static double mean(std::vector<float> const &vet);

    static float entropy_vector(std::vector<int> values);

    static float energy(std::vector<float> const &v);


};

#endif /* STATISTICS_H */
