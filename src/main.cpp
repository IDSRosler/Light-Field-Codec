#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <functional>
#include <memory>

#include "Typedef.h"
#include "utils.h"
#include "EncBitstreamWriter.h"
#include "LightField.h"
#include "LightFieldMetrics.h"
#include "Quantization.h"
#include "Timer.h"
#include "Transform.h"
#include "TextReport.h"
#include "EntropyEncoder.h"
#include "Prediction.h"


using namespace std;


vector<string> Split(const string &s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

float calculate_entropy(std::vector<int> values) {
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


void display_stage(const std::string& message) { std::cout << message << std::endl; }

int main(int argc, char **argv) {
    EncoderParameters encoderParameters;
    encoderParameters.parse_cli(argc, argv);
    encoderParameters.report();
    std::cout.flush();
    TextReport report;

#if STATISTICS_TIME
    Timer getBlock, rebuild, t, q, ti, qi, total_time;
#endif
#if STATISTICS_LOCAL
    Statistics statistics(encoderParameters.getPathOutput() + "statistics.csv");
    vector<float> input;
    int hypercube = 0;
#endif

    if (encoderParameters.display_stages)
        display_stage("[Loading light field]");

    LightField lf(encoderParameters.dim_LF, encoderParameters.getPathInput(),
                  encoderParameters.isLytro());

    const std::size_t SIZE = encoderParameters.dim_block.getNSamples();

    float orig4D[SIZE];
    float ti4D[SIZE];
    float qf4D[SIZE];
    float tf4D[SIZE];
    float qi4D[SIZE];
    float pf4D[SIZE];
    float pi4D[SIZE];
    int temp_lre[SIZE];

    static const std::vector<std::string> ch_names = {"Y", "Cb", "Cr"};


    Prediction predictor;
    //EDUARDO BEGIN
    Prediction newPredictor[3]{{(uint) ceil(encoderParameters.dim_LF.x / encoderParameters.dim_block.x)},
                               {(uint) ceil(encoderParameters.dim_LF.x / encoderParameters.dim_block.x)},
                               {(uint) ceil(encoderParameters.dim_LF.x / encoderParameters.dim_block.x)}};
    float ref4D[encoderParameters.dim_block.getNSamples()],
            res4D[encoderParameters.dim_block.getNSamples()];

    Transform transform(encoderParameters);

    LRE lre(encoderParameters.dim_block);
#if DPCM_DC
    DpcmDC dpcmDc[3]{{encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x}};
#endif

#if ENTROPY_TYPE
    EntropyEncoder encoder(&encoderParameters, 10000000);
#else
    EncBitstreamWriter encoder(&encoderParameters, 10000000);
#endif

    std::string sep = ",";

    const Point4D dimLF = encoderParameters.dim_LF;
    Point4D it_pos, dimBlock, stride_lf, stride_block;

    std::vector<index_t> seq_order(SIZE);
    std::iota(seq_order.begin(), seq_order.end(), 0);

#if STATISTICS_TIME
    total_time.tic();
#endif

    std::string tree;
    double ssim[3], psnr[3];
    double ssim_sum[3] = {0};
    double psnr_sum[3] = {0};
    double mean_ssim;
    double mean_psnr;

    volatile float total_steps = 3;
    total_steps *= std::ceil(dimLF.x / (float) encoderParameters.dim_block.x);
    total_steps *= std::ceil(dimLF.y / (float) encoderParameters.dim_block.y);
    total_steps *= std::ceil(dimLF.u / (float) encoderParameters.dim_block.u);
    total_steps *= std::ceil(dimLF.v / (float) encoderParameters.dim_block.v);
    float current_step = 1;

    if (encoderParameters.display_stages)
        display_stage("[Start encoding]");

    std::ofstream transform_stats;
    std::vector<index_t> scan_order;
    Point4D stride = make_stride(encoderParameters.dim_block);
    if (encoderParameters.export_transform_stats) {
        transform_stats = std::ofstream(encoderParameters.getPathOutput() + "transform_stats.csv");
        transform_stats << "X,Y,U,V,Channel,TransformDescriptor,SSE,"
                        << "ZerosBefore,EnergyBefore,EntropyBefore,"
                        << "ZerosAfter,EnergyAfter,EntropyAfter\n";
    }

#if STATISTICS_TIME
    total_time.tic();
#endif
    int block = 0;

    report.header({
                          "Position",
                          "Ch",
                          "Descriptor",
                          "RD-Cost",
                          "SSE",
                  });

    std::string transform_descriptor;
    double rd_cost;
    // angular
    for (it_pos.v = 0; it_pos.v < dimLF.v; it_pos.v += dimBlock.v) {
        for (it_pos.u = 0; it_pos.u < dimLF.u; it_pos.u += dimBlock.u) {
            // spatial
            for (it_pos.y = 0; it_pos.y < dimLF.y; it_pos.y += dimBlock.y) {
                for (it_pos.x = 0; it_pos.x < dimLF.x; it_pos.x += dimBlock.x) {

                    dimBlock = Point4D(
                            std::min(encoderParameters.dim_block.x, dimLF.x - it_pos.x),
                            std::min(encoderParameters.dim_block.y, dimLF.y - it_pos.y),
                            std::min(encoderParameters.dim_block.u, dimLF.u - it_pos.u),
                            std::min(encoderParameters.dim_block.v, dimLF.v - it_pos.v));

                    stride_lf =
                            Point4D(1, dimLF.x - dimBlock.x, dimLF.x * (dimLF.y - dimBlock.y),
                                    dimLF.x * dimLF.y * (dimLF.u - dimBlock.u));

                    stride_block = Point4D(1, encoderParameters.dim_block.x - dimBlock.x,
                                           (encoderParameters.dim_block.y - dimBlock.y) *
                                           encoderParameters.dim_block.x,
                                           (encoderParameters.dim_block.u - dimBlock.u) *
                                           encoderParameters.dim_block.x *
                                           encoderParameters.dim_block.y);


                    block++;

                    for (int it_channel = 0; it_channel < 3; ++it_channel) {
#if STATISTICS_TIME
                        getBlock.tic();
#endif
                        std::fill_n(orig4D, SIZE, 0);
                        std::fill_n(tf4D, SIZE, 0);
                        std::fill_n(qf4D, SIZE, 0);
                        std::fill_n(temp_lre, SIZE, 0);

                        lf.getBlock(orig4D, it_pos, dimBlock, stride_block,
                                    encoderParameters.dim_block, stride_lf, it_channel);

#if STATISTICS_TIME
                        getBlock.toc();
#endif // STATISTICS_TIME


#if STATISTICS_TIME
                        t.tic();
#endif // STATISTICS_TIME
                        if(encoderParameters.getPrediction() == "none"){
                            std::copy(orig4D, orig4D + SIZE, res4D);
                        } else if(encoderParameters.getPrediction() == "angular"){
                            newPredictor[it_channel].angularPredictionVector(it_pos.x, it_pos.y, orig4D,
                                                                       encoderParameters.dim_block, pf4D, block, ref4D);
                            //newPredictor[0].writeHeatMap(encoderParameters.getPathOutput());
                            newPredictor[it_channel].residuePred(orig4D, pf4D, encoderParameters.dim_block, res4D);
#if STATISTICS_LOCAL
                            input.clear();
                            for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                                input.push_back(res4D[i]);
                            }
                            statistics.compute_prediction_statistics(input);
                            statistics.write_prediction_statistics(hypercube, it_pos, dimBlock, ch_names[it_channel]);
#endif
                        } else if(encoderParameters.getPrediction() == "all"){
                            std::copy(orig4D, orig4D + SIZE, res4D);
                        } else if(encoderParameters.getPrediction() == "DC"){
                            std::copy(orig4D, orig4D + SIZE, res4D);
                        } else if(encoderParameters.getPrediction() == "IBC"){
                            std::copy(orig4D, orig4D + SIZE, res4D);
                        } else if(encoderParameters.getPrediction() == "SHIFT"){
                            std::transform(orig4D, orig4D + SIZE, res4D,
                                           [](auto value) { return value - 512; });
                        } else {
                            encoderParameters.prediction = "none";
                            std::copy(orig4D, orig4D + SIZE, res4D);
                        }

#if STATISTICS_TIME
                        t.toc();
#endif // STATISTICS_TIME
                        if (encoderParameters.enable_transforms) {
                            transform.set_position(it_channel, it_pos);
                            auto[_desc, _rd_cost] = transform.forward(res4D, qf4D, dimBlock);
                            transform_descriptor = _desc;
                            rd_cost = _rd_cost;
                        } else {
                            std::copy(res4D, res4D + SIZE, qf4D);
                        }
#if STATISTICS_TIME
                        t.tic();
#endif // STATISTICS_TIME

#if STATISTICS_TIME
                        t.toc();
#endif // STATISTICS_TIME

                        std::transform(qf4D, qf4D + SIZE, temp_lre,
                                       [](auto value) { return static_cast<int>(std::round(value)); });

                        auto lre_result = lre.encodeCZI(temp_lre, 0, SIZE);

#if ENTROPY_TYPE
                        encoder.encodeHypercube(temp_lre, encoderParameters.dim_block);
#else
                        auto lre_size = encoder.write4DBlock(temp_lre, SIZE, lre_result);
#endif

                        std::copy(temp_lre, temp_lre + SIZE, qi4D);


#if STATISTICS_TIME
                        ti.tic();
#endif
                        if (encoderParameters.enable_transforms) {
                            transform.inverse(transform_descriptor, qi4D, ti4D, dimBlock);
                        } else {
                            std::copy(qi4D, qi4D + SIZE, ti4D);
                        }

#if STATISTICS_TIME
                        ti.toc();
#endif // STATISTICS_TIME


                        if (encoderParameters.getPrediction() == "SHIFT") {
                            std::transform(ti4D, ti4D + SIZE, pi4D,
                                           [](auto value) { return value + 512; });
                        } else if(encoderParameters.getPrediction() != "none"){
                            newPredictor->recResiduePred(ti4D, pf4D, encoderParameters.dim_block, pi4D);
                        } else{
                            std::copy(ti4D, ti4D + SIZE, pi4D);
                        }


#if STATISTICS_TIME
                        rebuild.tic();
#endif // STATISTICS_TIME
                            lf.rebuild(pi4D, it_pos, dimBlock, stride_block, encoderParameters.dim_block, stride_lf,
                                       it_channel);

#if STATISTICS_TIME
                        rebuild.toc();
#endif // STATISTICS_TIME

                        //EDUARDO BEGIN

                        if(encoderParameters.getPrediction() != "none") {
                            newPredictor[it_channel].update(pi4D, true, encoderParameters.dim_block.getNSamples());
                        }
                        //EDUARDO END

                        encoder.write_completedBytes();


#if LFCODEC_EXPORT_MICROIMAGES_TRANSFORM
                         save_microimage(encoderParameters.getPathOutput(),
                                         it_pos,
                                         it_channel,
                                         qf4D,
                                         dimBlock,
                                         make_stride(Point4D(15, 15, 13, 13)),
                                         transform_descriptor);
#endif
                        if (encoderParameters.show_progress_bar)
                            progress_bar(current_step / total_steps, 50);

                        current_step++;
                        auto sse = std::inner_product(orig4D, orig4D + SIZE, pi4D, 0.0,
                                                      [](auto acc, auto val) { return acc + val; },
                                                      [](auto blk, auto rec) {
                                                          auto error = blk - rec;
                                                          return error * error;
                                                      });
                        if (encoderParameters.export_transform_stats) {
                            scan_order = generate_scan_order(dimBlock, stride);
                            // (x, y, u, v, ch, desc, rd_cost,
                            // Before#0, BeforEnergy, BeforEntropy,
                            // After#0, AfterEnergy, AfterEntropy,

                            int zeros_before = 0;
                            int zeros_after = 0;
                            float energy_before = 0;
                            float energy_after = 0;

                            std::vector<int> values_before;
                            std::vector<int> values_after;

                            for (auto i: scan_order) {
                                if (static_cast<int>(res4D[i]) == 0)  zeros_before++;
                                if (temp_lre[i] == 0) zeros_after++;
                                energy_before += res4D[i] * res4D[i];
                                energy_after += temp_lre[i] * temp_lre[i];
                                values_before.push_back(static_cast<int>(res4D[i]));
                                values_after.push_back(temp_lre[i]);
                            }
                            energy_before /= scan_order.size();
                            energy_after /= scan_order.size();

                            float entropy_before = calculate_entropy(values_before);
                            float entropy_after = calculate_entropy(values_after);

                            transform_stats << it_pos.x << sep
                                            << it_pos.y << sep
                                            << it_pos.u << sep
                                            << it_pos.v << sep
                                            << it_channel << sep
                                            << transform_descriptor << sep
                                            << sse << sep
                                            << zeros_before << sep
                                            << energy_before << sep
                                            << entropy_before << sep
                                            << zeros_after << sep
                                            << energy_after << sep
                                            << entropy_after << "\n";
                        }
                        if (encoderParameters.verbose) {
                            report.set_key("Position", it_pos);
                            report.set_key("Ch", ch_names[it_channel]);
                            report.set_key("Descriptor", transform_descriptor + transform.m_encoding_type);
                            report.set_key("SSE", sse);
                            report.set_key("RD-Cost", rd_cost);
                            report.print();
                        }

                        if (encoderParameters.export_blocks) {
                            const auto &path = encoderParameters.getPathOutput();

                            save_microimage(path, it_pos, it_channel, orig4D, dimBlock, stride, "_1_orig4D", 1);
                            save_microimage(path, it_pos, it_channel, res4D, dimBlock, stride, "_2_res4D", 1);
                            save_microimage(path, it_pos, it_channel, qf4D, dimBlock, stride, "_3_qf4D", 1);
                            save_microimage(path, it_pos, it_channel, qf4D, dimBlock, stride, "_3_qf4D_energy", 2);
                            save_microimage(path, it_pos, it_channel, qf4D, dimBlock, stride, "_3_qf4D_bitlength", 3);
                            save_microimage(path, it_pos, it_channel, ti4D, dimBlock, stride, "_4_ti4D", 1);
                            save_microimage(path, it_pos, it_channel, pf4D, dimBlock, stride, "_5_pf4D", 1);
                            save_microimage(path, it_pos, it_channel, pi4D, dimBlock, stride, "_6_pi4D", 1);
                        }
                    }
                    ++hypercube;
                }
            }
        }
    }

    if (encoderParameters.calculate_metrics) {
        if (encoderParameters.display_stages)
            display_stage("[Calculating quality metrics]");
        for (std::size_t v = 0; v < dimLF.v; v++)
            for (std::size_t u = 0; u < dimLF.u; u++)
                for (std::size_t ch = 0; ch < 3; ch++) {
                    auto ref = convertViewToMat(lf.yCbCr_original[ch], dimLF, u, v);
                    auto rec = convertViewToMat(lf.yCbCr[ch], dimLF, u, v);
                    psnr_sum[ch] += getPSNR(ref, rec, 10);
                    ssim_sum[ch] += getMSSIM(ref, rec, 10);
                }

        for (int ch = 0; ch < 3; ch++) {
            ssim[ch] = ssim_sum[ch] / static_cast<double>(dimLF.v * dimLF.u);
            psnr[ch] = psnr_sum[ch] / static_cast<double>(dimLF.v * dimLF.u);
        }

        mean_ssim = (6 * ssim[0] + ssim[1] + ssim[2]) / 8;
        mean_psnr = (6 * psnr[0] + psnr[1] + psnr[2]) / 8;
    }
#pragma clang optimize off
    if (encoderParameters.display_stages)
        display_stage("[Writing reconstructed Light Fields on disk]");
    lf.write(encoderParameters.getPathOutput());
#pragma clang optimize on
    encoder.finish_and_write();


    if (encoderParameters.display_stages)
        display_stage("[Done]");
#if STATISTICS_TIME
    total_time.toc();
#endif
    auto size_bytes = encoder.getTotalBytes();
    auto total_bpp = static_cast<double>(size_bytes)/static_cast<double>(dimLF.getNSamples());

    display_report(std::cout, "Total Bytes", size_bytes);
    display_report(std::cout, "Human readable size", to_human_readable(size_bytes));
    display_report(std::cout, "Bpp", total_bpp);
    if (encoderParameters.calculate_metrics) {
        display_report(std::cout, "PSNR-Y", psnr[0]);
        display_report(std::cout, "PSNR-U", psnr[1]);
        display_report(std::cout, "PSNR-V", psnr[2]);
        display_report(std::cout, "PSNR-YUV", mean_psnr);
        display_report(std::cout, "SSIM-Y", ssim[0]);
        display_report(std::cout, "SSIM-U", ssim[1]);
        display_report(std::cout, "SSIM-V", ssim[2]);
        display_report(std::cout, "SSIM-YUV", mean_ssim);
#if STATISTICS_TIME
        display_report(std::cout, "Transform::forward time: ", t.getTotalTime());
        display_report(std::cout, "Transform::inverse time: ", ti.getTotalTime());
#endif
    }
    cout << std::endl;

#if STATISTICS_LOCAL
    statistics.~Statistics();
#endif

#if STATISTICS_TIME
//  std::ofstream file_time;
//  file_time.open(encoderParameters.getPathOutput() + "time.csv");
//  file_time << "GetBlock Time(sec)" << sep << "T Time(sec)" << sep
//            << "Q Time(sec)" << sep << "TI Time(sec)" << sep << "QI Time(sec)"
//            << sep << "Rebuild Time(sec)" << sep << "Total Time(sec)" << sep;
//  file_time << std::endl;
//  file_time << getBlock.getTotalTime() << sep << t.getTotalTime() << sep
//            << q.getTotalTime() << sep << ti.getTotalTime() << sep
//            << qi.getTotalTime() << sep << rebuild.getTotalTime() << sep
//            << total_time.getTotalTime() << sep;
//  file_time << std::endl;
#endif

    return 0;
}
