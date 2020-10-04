#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <deque>
#include <algorithm>
#include <numeric>
#include <functional>

#include "DpcmDC.h"
#include "EncBitstreamWriter.h"
#include "EncoderParameters.h"
#include "LRE.h"
#include "LightField.h"
#include "LightFieldMetrics.h"
#include "Prediction.h"
#include "Quantization.h"
#include "ScanOrder.h"
#include "Statistics.h"
#include "Timer.h"
#include "Transform.h"
#include "Typedef.h"
#include "utils.h"
#include "TextReport.h"
using namespace std;


void display_stage(std::string message) { std::cout << message << std::endl; }

int main(int argc, char **argv) {
    EncoderParameters encoderParameters;
    encoderParameters.parse_cli(argc, argv);
    encoderParameters.report();
    std::cout.flush();
    TextReport report;

#if STATISTICS_TIME
    Timer getBlock, rebuild, t, q, ti, qi, total_time;
#endif

    if (encoderParameters.display_stages)
        display_stage("[Loading light field]");

    LightField lf(encoderParameters.dim_LF, encoderParameters.getPathInput(),
                  encoderParameters.isLytro());

    const std::size_t SIZE = encoderParameters.dim_block.getNSamples();

    std::shared_ptr<float[]> _orig4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _ti4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _qf4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _tf4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _qi4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _pf4D = std::shared_ptr<float[]>(new float[SIZE]);
    std::shared_ptr<float[]> _pi4D = std::shared_ptr<float[]>(new float[SIZE]);

    std::shared_ptr<int[]> _temp_lre = std::shared_ptr<int[]>(new int[SIZE]);

    auto orig4D = _orig4D.get();
    auto ti4D = _ti4D.get();
    auto qf4D = _qf4D.get();
    auto tf4D = _tf4D.get();
    auto qi4D = _qi4D.get();
    auto pf4D = _pf4D.get();
    auto pi4D = _pi4D.get();
    auto temp_lre = _temp_lre.get();

    static const std::vector<std::string> ch_names = {"Y", "Cb", "Cr"};

#if LFCODEC_USE_PREDICTION
    Prediction predictor;
#endif
    bool should_show_block = false;
    Transform transform(encoderParameters);

    // TODO: Entropy entropy(...)
    LRE lre(encoderParameters.dim_block);
#if DPCM_DC
    DpcmDC dpcmDc[3]{{encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x}};
#endif

    EncBitstreamWriter encoder(&encoderParameters, 10'000'000);
    std::string sep = ",";

    const Point4D dimLF = encoderParameters.dim_LF;
    Point4D it_pos, dimBlock, stride_lf, stride_block;

    Point4D blk_stride = make_stride(encoderParameters.dim_block);

    std::vector<index_t> seq_order(SIZE);
    std::iota(seq_order.begin(), seq_order.end(), 0);
    // for (std::size_t i = 0; i < SIZE; i++)
    //   seq_order[i] = i;

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

    std::map<char, int> transform_usage;

    auto z_order = generate_z_order_curve(encoderParameters.dim_block, blk_stride);

    report.header({
        "Position",
        "Ch",
        "TimerSubstr",
        "TimerFastFwd",
        "TimerFastInv",
        "TimerSplit",
        "TimerUniquePtr",
        "TimerStackCopy",
        "TimerRDCost",
        "TimerCacheCopy",
        "Descriptor"
    });

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
#endif
                        if ((it_pos.x / 15) == 41) {
                            // Dummy variable for breakpoint and debugging
                            volatile int xxx = 0;
                            // return 0;
                        }

#if LFCODEC_USE_PREDICTION
                        predictor.predict(orig4D, encoderParameters.dim_block, pf4D);
#else
                        std::copy(orig4D, orig4D + SIZE, pf4D);
                        for (std::size_t i = 0; i < SIZE; ++i) {
                          pf4D[i] = orig4D[i];
                        }
#endif

#if TRANSF_QUANT
                        transform.set_position(it_channel, it_pos);
#if STATISTICS_TIME
                        t.tic();
#endif

                        auto[descriptor, rd_cost] = transform.forward(pf4D, qf4D, dimBlock);


#if STATISTICS_TIME
                        t.toc();
#endif

                        std::transform(qf4D, qf4D + SIZE, temp_lre,
                                       [](auto value) { return static_cast<int>(value); });


                        auto lre_result = lre.encodeCZI(temp_lre, 0, SIZE);
                        auto lre_size = encoder.write4DBlock(temp_lre, SIZE, lre_result);



                        std::copy(temp_lre, temp_lre + SIZE, qi4D);

#else // TRANSF_QUANT
                        std::copy(pf4D, pf4D + SIZE, qf4D);
#endif // TRANSF_QUANT

#if TRANSF_QUANT

#if STATISTICS_TIME
                        ti.tic();
#endif

                        transform.inverse(descriptor, qi4D, ti4D, dimBlock);

#if STATISTICS_TIME
                        ti.toc();
#endif
#else // TRANSF_QUANT
                        std::copy(qi4D, qi4D + SIZE, ti4D);

#endif


#if LFCODEC_USE_PREDICTION
                        predictor.rec(ti4D, pi4D, dimBlock);
#else
                        std::copy(ti4D, ti4D + SIZE, pi4D);

#endif

#if STATISTICS_TIME
                        rebuild.tic();
#endif
                        lf.rebuild(pi4D, it_pos, dimBlock, stride_block,
                                   encoderParameters.dim_block, stride_lf, it_channel);

#if STATISTICS_TIME
                        rebuild.toc();
#endif

                        encoder.write_completedBytes();

                        if (should_show_block) {
                            show_block(it_channel, qf4D, dimBlock,
                                       make_stride(Point4D(15, 15, 13, 13)), "ref");
                        }
                        // save_microimage(encoderParameters.getPathOutput(),
                        //                 it_pos,
                        //                 it_channel,
                        //                 qf4D,
                        //                 dimBlock,
                        //                 make_stride(Point4D(15, 15, 13, 13)),
                        //                 descriptor);
                        if (encoderParameters.show_progress_bar)
                            progress_bar(current_step / total_steps, 50);

                        current_step++;

                        if (encoderParameters.verbose) {
                            report.set_key("Position", it_pos);
                            report.set_key("TimerSubstr", transform.m_timer_substr.getTotalTime() / transform.m_timer_substr.get_nCalls());
                            report.set_key("TimerFastFwd", transform.m_timer_md_forward_fast.getTotalTime() / transform.m_timer_md_forward_fast.get_nCalls());
                            report.set_key("TimerFastInv", transform.m_timer_md_inverse_fast.getTotalTime() / transform.m_timer_md_inverse_fast.get_nCalls());
                            report.set_key("TimerSplit", transform.m_timer_split.getTotalTime() / transform.m_timer_split.get_nCalls());
                            report.set_key("TimerUniquePtr", transform.m_timer_unique_ptr_alloc.getTotalTime() / transform.m_timer_unique_ptr_alloc.get_nCalls());
                            report.set_key("TimerStackCopy", transform.m_timer_stack_copy.getTotalTime() / transform.m_timer_stack_copy.get_nCalls());
                            report.set_key("TimerRDCost", transform.m_timer_rd_cost.getTotalTime() / transform.m_timer_rd_cost.get_nCalls());
                            report.set_key("TimerCacheCopy", transform.m_timer_cache_copy_overhead.getTotalTime() / transform.m_timer_cache_copy_overhead.get_nCalls());
                            report.set_key("Ch", ch_names[it_channel]);
                            report.set_key("Descriptor", descriptor + transform.m_encoding_type);
                            report.print();
                        }
                    }
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
            ssim[ch] = ssim_sum[ch] / (dimLF.v * dimLF.u);
            psnr[ch] = psnr_sum[ch] / (dimLF.v * dimLF.u);
        }

        mean_ssim = (6 * ssim[0] + ssim[1] + ssim[2]) / 8;
        mean_psnr = (6 * psnr[0] + psnr[1] + psnr[2]) / 8;
    }

    if (encoderParameters.display_stages)
        display_stage("[Writing reconstructed Light Fields on disk]");
    lf.write(encoderParameters.getPathOutput());
    encoder.finish_and_write();


    if (encoderParameters.display_stages)
        display_stage("[Done]");
#if STATISTICS_TIME
    total_time.toc();
#endif
    auto size_bytes = encoder.getTotalBytes();
    auto total_bpp = (float) (size_bytes * 1.25) / (float) dimLF.getNSamples();

    display_report(std::cout, "Total Bytes", size_bytes);
    display_report(std::cout, "Human readable size",
                   to_human_readable(size_bytes));
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
