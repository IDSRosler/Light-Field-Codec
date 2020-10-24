

#include "Statistics.h"

#include <cmath>
#include <cassert>
Statistics::Statistics(const std::string &filename) {
    file_out.open(filename);
    assert(file_out.is_open());
    //write_headers(file_out);
    this->write_headers();
} 

void Statistics::write_headers(){
    this->file_out
            << "hypercube"                  << sep
            << "channel"                    << sep
            << "pos_x"                      << sep
            << "pos_y"                      << sep
            << "pos_u"                      << sep
            << "pos_v"                      << sep
            << "bl_x"                       << sep
            << "bl_y"                       << sep
            << "bl_u"                       << sep
            << "bl_v"                       << sep
            << "entropy"                     << sep
            << "energy(x 0.00001)"                    << std::endl;
}

void Statistics::write_prediction_statistics(int hypercube, Point4D &pos, Point4D &dimBlock, std::string it_channel){
    this->file_out
            <<  hypercube                   << sep
            <<  it_channel                  << sep
            <<  pos.x                       << sep
            <<  pos.y                       << sep
            <<  pos.u                       << sep
            <<  pos.v                       << sep
            <<  dimBlock.x                  << sep
            <<  dimBlock.y                  << sep
            <<  dimBlock.u                  << sep
            <<  dimBlock.v                  << sep
            << this->v_entropy               << sep
            << this->v_energy * 0.00001              << std::endl;
}

void Statistics::write_headers(std::ofstream &output) {
    output 
        << "channel"                    << sep 
        << "pos_x"                      << sep
        << "pos_y"                      << sep
        << "pos_u"                      << sep
        << "pos_v"                      << sep
        << "bl_x"                       << sep
        << "bl_y"                       << sep
        << "bl_u"                       << sep
        << "bl_v"                       << sep
        << "segment"                    << sep
        << "dc"                         << sep           
        << "v_min"                      << sep 
        << "v_max"                      << sep 
        << "v_minAbs"                   << sep 
        << "v_maxAbs"                   << sep 
        << "zeros"                      << sep 
        << "ones"                       << sep 
        << "mean"                       << sep 
        << "median"                     << sep 
        << "std"                        << sep 
        << "variance"                   << sep 
        << "energy"                     << sep 
        << "entropy"                    << sep 
        << "sse"                        << sep       
        << "run_count"                  << sep 
        << "run_max"                    << sep 
        << "run_mean"                   << sep 
        << "run_std"                    << sep           
        << "run_var"                    << sep 
        << "posSO_last_nzero"           << sep 
        << "posSO_last_nzeroone"        << sep           
        << "pos_last_nzero"             << sep 
        << "pos_last_nzeroone"          << sep 
        << "bits_per_4D_Block"          << std::endl;
}

void Statistics::write(Point4D &pos, Point4D &dimBlock, std::size_t it_channel, std::string segment)
{
    write(file_out, pos, dimBlock, it_channel, segment);
}
void Statistics::write(std::ostream &output, Point4D &pos, Point4D &dimBlock, std::size_t it_channel, std::string segment) {
    std::stringstream ss;
    ss << "\"" << segment << "\"";
    output 
        << it_channel                   << sep 
        << pos.x                        << sep 
        << pos.y                        << sep 
        << pos.u                        << sep 
        << pos.v                        << sep
        << dimBlock.x                   << sep 
        << dimBlock.y                   << sep 
        << dimBlock.u                   << sep 
        << dimBlock.v                   << sep
        << ss.str()                     << sep
        << this->dc                     << sep 
        << this->v_min                  << sep 
        << this->v_max                  << sep 
        << this->v_minAbs               << sep
        << this->v_maxAbs               << sep 
        << this->num_zeros              << sep 
        << this->num_ones               << sep
        << this->v_mean                 << sep 
        << this->v_median               << sep 
        << this->v_std                  << sep 
        << this->v_variance             << sep 
        << this->v_energy               << sep 
        << this->v_entropy              << sep 
        << this->sse                    << sep
        << this->cov                    << sep 
        << std::endl;
}

void Statistics::write(Point4D &pos,
                       Point4D &dimBlock,
                       std::size_t it_channel,
                       std::vector<LRE_struct> &lre_result,
                       std::size_t bits_per_4D_Block) {
    std::ostream &output = (this->file_out.is_open()) ? this->file_out : std::cout;

    std::size_t max_run = 0;

    std::vector<float> run_zeros;

    for (auto it_result: lre_result) {
        max_run = (it_result.run > max_run) ? it_result.run : max_run;
        run_zeros.push_back(it_result.run);
    }

    this->run_mean = Statistics::mean(run_zeros);
    this->run_var = Statistics::variance(run_zeros, this->run_mean);
    this->run_std = sqrt(this->run_var);

    output 
        << it_channel                  << sep 
        << pos.x                       << sep 
        << pos.y                       << sep 
        << pos.u                       << sep 
        << pos.v                       << sep 
        << dimBlock.x                  << sep 
        << dimBlock.y                  << sep 
        << dimBlock.u                  << sep 
        << dimBlock.v                  << sep 
        << ""                          << sep
        << this->dc                    << sep 
        << this->v_min                 << sep 
        << this->v_max                 << sep 
        << this->v_minAbs              << sep        
        << this->v_maxAbs              << sep 
        << this->num_zeros             << sep 
        << this->num_ones              << sep 
        << this->v_mean                << sep 
        << this->v_median              << sep           
        << this->v_std                 << sep 
        << this->v_variance            << sep 
        << this->v_energy              << sep           
        << this->v_entropy             << sep 
        << this->sse                   << sep 
        << lre_result.size()           << sep 
        << max_run                     << sep 
        << this->run_mean              << sep 
        << this->run_std               << sep 
        << this->run_var               << sep 
        << this->posSO_last_nzero      << sep 
        << this->posSO_last_nzeroone   << sep 
        << this->pos_last_nzero        << sep 
        << this->pos_last_nzeroone     << sep 
        << bits_per_4D_Block           << std::endl;
}

Statistics::~Statistics() {
    if (this->file_out.is_open()) { this->file_out.close(); }
}

void Statistics::compute_prediction_statistics(const std::vector<float> &input){
    this->v_energy = Statistics::energy(input);
    this->v_entropy = Statistics::entropy_vector(std::vector<int>(input.begin(), input.end()));
}

void Statistics::compute(const std::vector<float> &input, const ushort *scan_order) {
    std::vector<float> values;
    this->dc = input[scan_order ? scan_order[0] : 0];
    this->num_zeros = this->num_ones = 0;
    this->posSO_last_nzero = this->posSO_last_nzeroone = -1;
    this->pos_last_nzero = this->pos_last_nzeroone = -1;

    for (std::size_t i = 0; i < input.size(); ++i) {
        float value = input[scan_order ? scan_order[i] : i];
        assert(!std::isnan(value));
        if (std::round(value) == 0) {
            ++this->num_zeros;
        } else if (std::round(value) == 1) { // epsilon = 0.1
            ++this->num_ones;
            this->posSO_last_nzero = scan_order ? scan_order[i] : i;
            this->pos_last_nzero = i;
        } else {
            values.push_back(value);
            this->posSO_last_nzero = this->posSO_last_nzeroone = scan_order ? scan_order[i] : i;
            this->pos_last_nzero = this->pos_last_nzeroone = i;
        }
    }

    this->v_mean = Statistics::mean(values);
    // BUG: values arent sorted
    this->v_median = Statistics::median(values);
    this->v_variance = Statistics::variance(values, this->v_mean);
    this->v_std = sqrt(this->v_variance);

    this->v_energy = Statistics::energy(input);
    this->v_entropy = Statistics::entropy_vector(std::vector<int>(input.begin(), input.end()));

    assert(!std::isnan(this->v_mean));

    this->v_min = this->v_max = this->v_minAbs = this->v_maxAbs = 0;
    if (!values.empty()) {
        this->v_min = Statistics::min_vet(std::vector<float>(values.begin() + 1, values.end()));
        this->v_max = Statistics::max_vet(std::vector<float>(values.begin() + 1, values.end()));
        this->v_minAbs =
        Statistics::minAbs_vet(std::vector<float>(values.begin() + 1, values.end()));
        this->v_maxAbs =
        Statistics::maxAbs_vet(std::vector<float>(values.begin() + 1, values.end()));
    }
}

float Statistics::compute_sse(float *orig,
                             float *ref,
                             const Point4D &dim_block,
                             const Point4D &stride_block) {
    float error;
    this->sse = 0;
    this->cov = 0;
    float *it_orig = orig, *it_ref = ref;

    int count = 0;
    for (std::size_t it_v = 0; it_v < dim_block.v; ++it_v) {
        for (std::size_t it_u = 0; it_u < dim_block.u; ++it_u) {
            for (std::size_t it_y = 0; it_y < dim_block.y; ++it_y) {
                for (std::size_t it_x = 0; it_x < dim_block.x; ++it_x) {
                    // std::size_t pos = it_x + it_y * 15 + it_u * 15 * 15 + it_v * 15 * 15 * 15;

                    error = *it_ref - *it_orig;
                    this->sse += error * error;
                    ++count;
                    it_orig += stride_block.x;
                    it_ref += stride_block.x;
                }
                it_orig += stride_block.y;
                it_ref += stride_block.y;
            }
            it_orig += stride_block.u;
            it_ref += stride_block.u;
        }
        it_orig += stride_block.v;
        it_ref += stride_block.v;
    }
    return sse;
}

float Statistics::entropy_vector(std::vector<int> values) {
    std::vector<std::tuple<int, std::size_t, float>> elements_count;

    std::sort(values.begin(), values.end());
    auto unique_values = values;

    unique_values.resize(
    std::distance(unique_values.begin(), std::unique(unique_values.begin(), unique_values.end())));

    for (auto &i: unique_values) {
        int v_count = std::count(values.begin(), values.end(), i);
        float v_freq = float(v_count) / values.size();
        elements_count.emplace_back(i, v_count, v_freq);
    }

    float prob, entropy = 0;
    for (auto &it: elements_count) {
        prob = std::get<2>(it);
        entropy += (float)prob * std::log2(prob);
    }
    return -entropy;
}

float Statistics::energy(std::vector<float> const &v) {
    float power = 0.0f;
    for (auto &i: v) {
        power += (float)i * (float)i;
    }
    return power;
}

double Statistics::variance(std::vector<float> const &vet, double mean) {
    double variance = 0;

    for (auto it: vet) {
        variance += pow((it - mean), 2);
    }

    return (variance / vet.size());
}

double Statistics::median(std::vector<float> vet) {
    size_t size = vet.size();

    if (size == 0) return 0;

    return (size % 2 == 0) ? ((vet[size / 2 - 1] + vet[size / 2]) / 2) : (vet[size / 2]);
}

double Statistics::mean(std::vector<float> const &vet) {
    if (vet.empty()) return 0;

    float sum = 0.0;
    float old_sum = 0.0;
    for (auto &i: vet) {
        sum += i;
        assert(!std::isnan(sum));
        old_sum = sum;
    }

    return sum / (float)vet.size();
}

double Statistics::get_mean() const { return this->v_mean; }

double Statistics::get_median() const { return this->v_median; }

double Statistics::get_variance() const { return this->v_variance; }

double Statistics::get_std() const { return this->v_std; }

unsigned long Statistics::get_zeros() const { return this->num_zeros; }

unsigned long Statistics::get_ones() const { return this->num_ones; }

double Statistics::get_entropy() const { return this->v_entropy; }

double Statistics::get_power() const { return this->v_energy; }

float Statistics::min_vet(const std::vector<float> &input) {
    if (input.empty()) return 0;
    float t_min = input[0];
    for (auto it: input)
        t_min = std::min(it, t_min);
    return t_min;
}

float Statistics::max_vet(const std::vector<float> &input) {
    if (input.empty()) return 0;
    float t_max = input[0];
    for (auto it: input)
        t_max = std::max(it, t_max);
    return t_max;
}

float Statistics::minAbs_vet(const std::vector<float> &input) {
    if (input.empty()) return 0;
    auto t_min_abs = std::abs(input[0]);
    for (auto it: input)
        t_min_abs = std::min(std::abs(it), t_min_abs);
    return t_min_abs;
}

float Statistics::maxAbs_vet(const std::vector<float> &input) {
    if (input.empty()) return 0;
    auto t_max_abs = std::abs(input[0]);
    for (auto it: input)
        t_max_abs = std::max(std::abs(it), t_max_abs);
    return t_max_abs;
}

// std::map<float, int> Statistics::calculate_pdf(const std::vector<float> &input, unsigned range) {
//     std::map<float, int> hist;
//     for (const auto &val: input) {
//             auto key = std::trunc(val);
//             if (range > 0)
//                 key /= range;
//             auto old = hist[key];
//             hist[key] = old + 1;
//         }
//     return hist;
// }