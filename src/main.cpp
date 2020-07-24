#include <iostream>
#include <stdexcept>
#include <string>
#include <map>

#include "DpcmDC.h"
#include "EncBitstreamWriter.h"
#include "EncoderParameters.h"
#include "LRE.h"
#include "LightField.h"
#include "Prediction.h"
#include "Quantization.h"
#include "Statistics.h"
#include "Time.h"
#include "Transform.h"
#include "Typedef.h"
#include "utils.h"

using namespace std;

int main(int argc, char **argv) {
    EncoderParameters encoderParameters;
    encoderParameters.parse_cli(argc, argv);
    encoderParameters.report();

#if STATISTICS_TIME
    Time getBlock, rebuild, t, q, ti, qi, total_time;
#endif

    LightField lf(encoderParameters.dim_LF, encoderParameters.getPathInput(),
                  encoderParameters.isLytro());

    float orig4D[encoderParameters.dim_block.getNSamples()],
    ti4D[encoderParameters.dim_block.getNSamples()],
    qf4D[encoderParameters.dim_block.getNSamples()],
    tf4D[encoderParameters.dim_block.getNSamples()],
    qi4D[encoderParameters.dim_block.getNSamples()],
    pf4D[encoderParameters.dim_block.getNSamples()],
    pi4D[encoderParameters.dim_block.getNSamples()];

    int temp_lre[encoderParameters.dim_block.getNSamples()];
    uint bits_per_4D_Block = 0;

#if LFCODEC_USE_PREDICTION
    Prediction predictor;
#endif
    volatile bool should_show_block = false;
    Transform transform(encoderParameters);
    auto transform_type = Transform::get_type(encoderParameters.transform);
    
    // TODO: Entropy entropy(...)
    LRE lre(encoderParameters.dim_block.getNSamples() == 15 * 15 * 15 * 15);
    // todo: fixed block size (15x15x15x15) or (15x15x13x13)
#if DPCM_DC
    DpcmDC dpcmDc[3] {{encoderParameters.dim_LF.x},
                      {encoderParameters.dim_LF.x},
                      {encoderParameters.dim_LF.x}};
#endif

    EncBitstreamWriter encoder(&encoderParameters, 10'000'000);
    // encoder.writeHeader();

#if STATISTICS_LOCAL
    
    Statistics statistics_qf(encoderParameters.getPathOutput() +
                             "localStatistics_quantization.csv");
#endif

    std::string sep = ",";

#if TRACE_TRANSF
    std::ofstream file_traceTransf;
    file_traceTransf.open(encoderParameters.getPathOutput() + "traceTransf.csv");
    file_traceTransf << "channel" << sep << "pos_x" << sep << "pos_y" << sep << "pos_u" << sep
                     << "pos_v" << sep << "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v"
                     << sep;

    for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
        file_traceTransf << i;
        if (i != encoderParameters.dim_block.getNSamples() - 1) { file_traceTransf << sep; }
    }

    file_traceTransf << std::endl;
#endif

#if TRACE_QUANT
    std::ofstream file_traceQuant;
    file_traceQuant.open(encoderParameters.getPathOutput() + "traceQuant.csv");
    file_traceQuant << "channel" << sep << "pos_x" << sep << "pos_y" << sep << "pos_u" << sep
                    << "pos_v" << sep << "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v"
                    << sep;

    for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
        file_traceQuant << i;
        if (i != encoderParameters.dim_block.getNSamples() - 1) { file_traceQuant << sep; }
    }

    file_traceQuant << std::endl;
#endif

#if TRACE_LRE
    std::ofstream file_traceLRE;
    file_traceLRE.open(encoderParameters.getPathOutput() + "traceLre.csv");
    file_traceLRE << "channel" << sep << "pos_x" << sep << "pos_y" << sep << "pos_u" << sep
                  << "pos_v" << sep << "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v"
                  << sep << "level" << sep << "run" << std::endl;
#endif
    const Point4D dimLF = encoderParameters.dim_LF;

    Point4D it_pos, dimBlock, stride_lf, stride_block;

#if LFCODEC_TRANSFORM_HISTOGRAM
    std::map<float, int> histogram;
    Point4D stride = make_stride(encoderParameters.dim_block);
    std::ofstream file_histogram;
    file_histogram.open(encoderParameters.getPathOutput() + "histogram.csv");
    file_histogram << "value" << sep << "frequency" << std::endl;
#endif

#if STATISTICS_TIME
    total_time.tic();
#endif
    std::string tree;

    volatile float total_steps = 3;
    total_steps *= std::ceil(dimLF.x / (float)encoderParameters.dim_block.x);
    total_steps *= std::ceil(dimLF.y / (float)encoderParameters.dim_block.y);
    total_steps *= std::ceil(dimLF.u / (float)encoderParameters.dim_block.u);
    total_steps *= std::ceil(dimLF.v / (float)encoderParameters.dim_block.v);
    float current_step = 1;

    // angular
    for (it_pos.v = 0; it_pos.v < dimLF.v; it_pos.v += dimBlock.v) {
        for (it_pos.u = 0; it_pos.u < dimLF.u; it_pos.u += dimBlock.u) {
            // spatial
            for (it_pos.y = 0; it_pos.y < dimLF.y; it_pos.y += dimBlock.y) {
                for (it_pos.x = 0; it_pos.x < dimLF.x; it_pos.x += dimBlock.x) {
                    dimBlock = Point4D(std::min(encoderParameters.dim_block.x, dimLF.x - it_pos.x),
                                       std::min(encoderParameters.dim_block.y, dimLF.y - it_pos.y),
                                       std::min(encoderParameters.dim_block.u, dimLF.u - it_pos.u),
                                       std::min(encoderParameters.dim_block.v, dimLF.v - it_pos.v));

                    stride_lf = Point4D(1, dimLF.x - dimBlock.x, dimLF.x * (dimLF.y - dimBlock.y),
                                        dimLF.x * dimLF.y * (dimLF.u - dimBlock.u));

                    stride_block = Point4D(
                    1, encoderParameters.dim_block.x - dimBlock.x,
                    (encoderParameters.dim_block.y - dimBlock.y) * encoderParameters.dim_block.x,
                    (encoderParameters.dim_block.u - dimBlock.u) * encoderParameters.dim_block.x *
                    encoderParameters.dim_block.y);

                    // std::cout << "Pos: " << it_pos << "\tBlock Size: " << dimBlock << std::endl;

                    

                    for (int it_channel = 0; it_channel < 3; ++it_channel) {
#if STATISTICS_TIME
                        getBlock.tic();
#endif
                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i)
                            orig4D[i] = tf4D[i] = qf4D[i] = 0;

                        lf.getBlock(orig4D, it_pos, dimBlock, stride_block,
                                    encoderParameters.dim_block, stride_lf, it_channel);

#if STATISTICS_TIME
                        getBlock.toc();
#endif
                        if (it_pos.x + 15 > dimLF.x)
                        {
                            // Dummy variable for breakpoint and debugging
                            volatile int xxx = 0;
                        }
#if TRANSF_QUANT

#if STATISTICS_TIME
                        t.tic();
#endif
#if LFCODEC_USE_PREDICTION
                        predictor.predict(orig4D, encoderParameters.dim_block, pf4D);
#else
                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                            pf4D[i] = orig4D[i];
                        }
#endif
                        
                        transform.set_position(it_channel, it_pos);
                        tree = transform.forward(transform_type, pf4D, qf4D, dimBlock);

                        if (false) {

                        printf("Pos(x=%02d,y=%02d,u=%02d,v=%02d,ch=%d) tree=%s\n", 
                               it_pos.x / 15, it_pos.y / 15, it_pos.u / 15, it_pos.v / 15, it_channel,
                               tree.c_str());
                        }

#if STATISTICS_TIME
                        t.toc();
#endif

#if TRACE_TRANSF
                        file_traceTransf << it_channel << sep <<

                        it_pos.x << sep << it_pos.y << sep << it_pos.u << sep << it_pos.v << sep <<

                        dimBlock.x << sep << dimBlock.y << sep << dimBlock.u << sep << dimBlock.v
                                         << sep;

                        for (auto it: tf4D) {
                            file_traceTransf << it << sep;
                        }

                        file_traceTransf << std::endl;
#endif
#if LFCODEC_TRANSFORM_HISTOGRAM
                        for (int v = 0; v < dimBlock.v; v++)
                            for (int u = 0; u < dimBlock.u; u++)
                                for (int y = 0; y < dimBlock.y; y++)
                                    for (int x = 0; x < dimBlock.x; x++) {
                                        auto index = offset(x, y, u, v, stride);
                                        auto value = tf4D[index];
                                        histogram[std::trunc(value)] += 1;
                                    }
#endif
#if QUANTIZATION

#if STATISTICS_TIME
                        q.tic();
#endif

                        // quantization.foward(Quantization::HOMOGENEOUS, tf4D, qf4D);

#if DPCM_DC
                        qf4D[0] -= (float)dpcmDc[it_channel].get_reference(it_pos.x, it_pos.y);
#endif

#if STATISTICS_TIME
                        q.toc();
#endif
#if TRACE_QUANT
                        file_traceQuant << it_channel << sep << it_pos.x << sep << it_pos.y << sep
                                        << it_pos.u << sep << it_pos.v << sep <<

                        dimBlock.x << sep << dimBlock.y << sep << dimBlock.u << sep << dimBlock.v
                                        << sep;

                        for (auto it: qf4D) {
                            file_traceQuant << (int)it << sep;
                        }

                        file_traceQuant << std::endl;
#endif
#else  // QUANTIZATION
                        std::copy(tf4D, tf4D + encoderParameters.dim_block.getNSamples(), qf4D);
#endif // QUANTIZATION

                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                            temp_lre[i] = (int)std::trunc(qf4D[i]);
                        }

                        auto lre_result =
                        lre.encodeCZI(temp_lre, 0, encoderParameters.dim_block.getNSamples());
#if TRACE_LRE
                        for (auto it_lre: lre_result) {
                            file_traceLRE << it_channel << sep << it_pos.x << sep << it_pos.y << sep
                                          << it_pos.u << sep << it_pos.v << sep <<

                            dimBlock.x << sep << dimBlock.y << sep << dimBlock.u << sep
                                          << dimBlock.v << sep <<

                            it_lre.level << sep << it_lre.run << std::endl;
                        }
#endif
                        bits_per_4D_Block = encoder.write4DBlock(
                        temp_lre, encoderParameters.dim_block.getNSamples(), lre_result);

#if STATISTICS_TIME
                        qi.tic();
#endif
#if QUANTIZATION
#if DPCM_DC
                        qf4D[0] += (float)dpcmDc[it_channel].get_reference(it_pos.x, it_pos.y);
#endif
                        // quantization.inverse(Quantization::HOMOGENEOUS, qf4D, qi4D);

#if STATISTICS_TIME
                        qi.toc();
#endif

#else /* NO_QUANTIZATION */
                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                            qi4D[i] = temp_lre[i];
                            // qi4D[i] = qf4D[i];
                        }
#endif

#if STATISTICS_TIME
                        ti.tic();
#endif


                        transform.inverse(tree, qf4D, ti4D, dimBlock);
#if LFCODEC_USE_PREDICTION
                        predictor.rec(ti4D, pi4D, dimBlock);
#else
                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                            pi4D[i] = ti4D[i];
                        }
#endif

#if STATISTICS_TIME
                        ti.toc();
#endif

#else /* NO_TRANSF_QUANT */
                        std::copy(orig4D, orig4D + dimBlock.getNSamples(), ti4D);
#endif

#if STATISTICS_TIME
                        rebuild.tic();
#endif
                        // lf.rebuild(ti4D, it_pos, dimBlock, stride_block,
                        // encoderParameters.dim_block, stride_lf, it_channel);
                        lf.rebuild(pi4D, it_pos, dimBlock, stride_block,
                                   encoderParameters.dim_block, stride_lf, it_channel);

#if STATISTICS_TIME
                        rebuild.toc();
#endif

#if STATISTICS_LOCAL
#if TRANSF_QUANT
#if QUANTIZATION

                        statistics_qf.compute(
                        std::vector<float>(qf4D, qf4D + encoderParameters.dim_block.getNSamples()),
                        lre.getScanOrder());
                        statistics_qf.compute_sse(orig4D, qi4D, dimBlock, stride_block);
                        statistics_qf.write(it_pos, dimBlock, it_channel, lre_result,
                                            bits_per_4D_Block);
#endif // QUANTIZATION
#endif // TRANSF_QUANT
#endif // STATISTICS_LOCAL

#if DPCM_DC
                        dpcmDc[it_channel].update((int)qf4D[0], true);
#endif
                        encoder.write_completedBytes();

                        if (should_show_block) {
                            show_block(it_channel, orig4D, dimBlock, make_stride(Point4D(15,15,13,13)), "ref");
                            show_block(it_channel, pi4D, dimBlock, make_stride(Point4D(15,15,13,13)), "rec");
                        }
                        
                        if (true)
                            progress_bar(current_step++ / total_steps, 50);
                    }
                }
            }
        }
    }

    lf.write(encoderParameters.getPathOutput());
    encoder.finish_and_write();
    encoder.~EncBitstreamWriter();

#if STATISTICS_TIME
    total_time.toc();
#endif

    cout << "\n#########################################################" << endl;
    cout << "Total bytes:\t" << encoder.getTotalBytes() << endl;
    cout << "#########################################################" << endl;

#if STATISTICS_GLOBAL
// TODO: statistics (global)
#endif

#if TRACE_QUANT
    file_traceQuant.close();
#endif

#if TRACE_LRE
    file_traceLRE.close();
#endif
#if LFCODEC_TRANSFORM_HISTOGRAM
    for (const auto &[val, freq]: histogram)
        file_histogram << val << sep << freq << std::endl;
    file_histogram.close();
#endif
#if STATISTICS_TIME
    std::ofstream file_time;
    file_time.open(encoderParameters.getPathOutput() + "time.csv");
    file_time << "GetBlock Time(sec)" << sep << "T Time(sec)" << sep << "Q Time(sec)" << sep
              << "TI Time(sec)" << sep << "QI Time(sec)" << sep << "Rebuild Time(sec)" << sep
              << "Total Time(sec)" << sep;
    file_time << std::endl;
    file_time << getBlock.getTotalTime() << sep << t.getTotalTime() << sep << q.getTotalTime()
              << sep << ti.getTotalTime() << sep << qi.getTotalTime() << sep
              << rebuild.getTotalTime() << sep << total_time.getTotalTime() << sep;
    file_time << std::endl;
#endif

    return 0;
}