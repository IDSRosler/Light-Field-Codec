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
#include <map>

#include "Point4D.h"
#include "LRE.h"


class Statistics {
public:
    float epsilon = 1e-1;

    Statistics() = default;

    explicit Statistics(const std::string &file);

    virtual ~Statistics();

    void write_headers(std::ofstream &output);
    void write_headers();

    void write(Point4D &pos, Point4D &dimBlock, std::size_t it_channel, std::vector<LRE_struct> &lre_result, std::size_t bits_per_4D_Block);

    void write(Point4D &pos, Point4D &dimBlock, std::size_t it_channel, std::string segment = "");
    void write(std::ostream &output, Point4D &pos, Point4D &dimBlock, std::size_t it_channel, std::string segment = "");

    void write_prediction_statistics(int hypercube, Point4D &pos, Point4D &dimBlock, std::string it_channel);
    void compute_prediction_statistics(const std::vector<float> &input);

    double get_mean() const;

    double get_median() const;

    double get_variance() const;

    double get_std() const;

    unsigned long get_zeros() const;

    unsigned long get_ones() const;

    double get_entropy() const;

    double get_power() const;

    void compute(const std::vector<float> &input, const ushort *scan_order);

    float compute_sse(float *orig, float *ref, const Point4D &dim_block, const Point4D &stride_block);

    static float min_vet(const std::vector<float> &input);

    static float max_vet(const std::vector<float> &input);

    static float minAbs_vet(const std::vector<float> &input);

    static float maxAbs_vet(const std::vector<float> &input);

    // static std::map<float, int> calculate_pdf(const std::vector<float> &input, unsigned range);

private:
    std::string sep{","};

    std::ofstream file_out;
    unsigned long num_zeros{}, num_ones{};
    double v_mean{}, v_median{}, v_variance{}, v_std{};
    double v_entropy{}, v_energy{};
    float dc, sse{0}, cov{0}, v_max, v_min, v_maxAbs, v_minAbs;

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
