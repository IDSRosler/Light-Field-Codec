#include <iostream>

#include "EncoderParameters.h"
#include "Quantization.h"
#include "Transform.h"
#include "LightField.h"
#include "Statistics.h"
#include "LRE.h"
#include "DpcmDC.h"
#include "EncBitstreamWriter.h"
#include "Typedef.h"
#include "Time.h"
#include "Prediction.h"
#include "Tree.h"

using namespace std;

template<class T>
void printVector(const std::string &msg, T *vet, int size) {
    std::cout << std::endl << msg << std::endl;
    for (int i = 0; i < size; ++i) {
        std::cout << vet[i] << "\t";
    }
    std::cout << std::endl;
}

#if HEXADECA_TREE
vector<string> Split(const string& s, char delimiter) {
    vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
#endif

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

    //Prediction predictor;
    Transform transform(encoderParameters.dim_block);
    Quantization quantization(encoderParameters.dim_block, encoderParameters.getQp(),
                              encoderParameters.quant_weight_100);
    // TODO: Entropy entropy(...)
    LRE lre(encoderParameters.dim_block.getNSamples() ==
            15 * 15 * 15 * 15); //todo: fixed block size (15x15x15x15) or (15x15x13x13)
#if DPCM_DC
    DpcmDC dpcmDc[3]{{encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x},
                     {encoderParameters.dim_LF.x}};
#endif

    EncBitstreamWriter encoder(&encoderParameters, 10000000);
    //encoder.writeHeader();

#if STATISTICS_LOCAL
    //    Statistics statistics_tf(encoderParameters.getPathOutput() + "localStatistics_transform.csv");
        Statistics statistics_qf(encoderParameters.getPathOutput() + "localStatistics_quantization.csv");
#endif

    std::string sep = ",";

#if HEXADECA_TREE
    uint hypercubo = 0;
    Tree tree;
    Node *root = nullptr;
    string light_field_name = Split(encoderParameters.getPathInput(), '/').back();
    tree.OpenFile(encoderParameters.getPathOutput());
#endif

#if TRACE_TRANSF
    std::ofstream file_traceTransf;
    file_traceTransf.open(encoderParameters.getPathOutput() + "traceTransf.csv");
    file_traceTransf <<
                     "channel" << sep <<
                     "pos_x" << sep << "pos_y" << sep << "pos_u" << sep << "pos_v" << sep <<
                     "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v" << sep;

    for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
        file_traceTransf << i;
        if (i != encoderParameters.dim_block.getNSamples() - 1) {
            file_traceTransf << sep;
        }
    }
    file_traceTransf << std::endl;
#endif

#if TRACE_QUANT
    std::ofstream file_traceQuant;
    file_traceQuant.open(encoderParameters.getPathOutput() + "traceQuant.csv");
    file_traceQuant <<
                    "channel" << sep <<
                    "pos_x" << sep << "pos_y" << sep << "pos_u" << sep << "pos_v" << sep <<
                    "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v" << sep;

    for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
        file_traceQuant << i;
        if (i != encoderParameters.dim_block.getNSamples() - 1) {
            file_traceQuant << sep;
        }
    }
    file_traceQuant << std::endl;
#endif

#if TRACE_LRE
    std::ofstream file_traceLRE;
    file_traceLRE.open(encoderParameters.getPathOutput() + "traceLre.csv");
    file_traceLRE <<
                  "channel" << sep <<
                  "pos_x" << sep << "pos_y" << sep << "pos_u" << sep << "pos_v" << sep <<
                  "bl_x" << sep << "bl_y" << sep << "bl_u" << sep << "bl_v" << sep <<
                  "level" << sep << "run" << std::endl;
#endif
    const Point4D dimLF = encoderParameters.dim_LF;

    Point4D it_pos, dimBlock, stride_lf, stride_block;

#if STATISTICS_TIME
    total_time.tic();
#endif

#if HEXADECA_TREE
    Point_4D hypercubo_pos;

    for (it_pos.v = 0, hypercubo_pos.v = 0; it_pos.v < dimLF.v; it_pos.v += dimBlock.v, hypercubo_pos.v++) { // angular
        for (it_pos.u = 0, hypercubo_pos.u = 0; it_pos.u < dimLF.u; it_pos.u += dimBlock.u, hypercubo_pos.u++) {

            for (it_pos.y = 0, hypercubo_pos.y = 0; it_pos.y < dimLF.y; it_pos.y += dimBlock.y, hypercubo_pos.y++) { // spatial
                for (it_pos.x = 0, hypercubo_pos.x = 0; it_pos.x < dimLF.x; it_pos.x += dimBlock.x, hypercubo_pos.x++) {
#else

    for (it_pos.v = 0; it_pos.v < dimLF.v; it_pos.v += dimBlock.v) { // angular
        for (it_pos.u = 0; it_pos.u < dimLF.u; it_pos.u += dimBlock.u) {

            for (it_pos.y = 0; it_pos.y < dimLF.y; it_pos.y += dimBlock.y) { // spatial
                for (it_pos.x = 0; it_pos.x < dimLF.x; it_pos.x += dimBlock.x) {
#endif

                    dimBlock = Point4D(std::min(encoderParameters.dim_block.x, dimLF.x - it_pos.x),
                                       std::min(encoderParameters.dim_block.y, dimLF.y - it_pos.y),
                                       std::min(encoderParameters.dim_block.u, dimLF.u - it_pos.u),
                                       std::min(encoderParameters.dim_block.v, dimLF.v - it_pos.v)
                    );

                    stride_lf = Point4D(1, dimLF.x - dimBlock.x,
                                        dimLF.x * (dimLF.y - dimBlock.y),
                                        dimLF.x * dimLF.y * (dimLF.u - dimBlock.u)
                    );

                    stride_block = Point4D(1, encoderParameters.dim_block.x - dimBlock.x,
                                           (encoderParameters.dim_block.y - dimBlock.y) *
                                           encoderParameters.dim_block.x,
                                           (encoderParameters.dim_block.u - dimBlock.u) *
                                           encoderParameters.dim_block.x * encoderParameters.dim_block.y);

                    std::cout << "Pos: " << it_pos << "\tBlock Size: " << dimBlock << std::endl;
                    for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                        orig4D[i] = tf4D[i] = qf4D[i] = 0;
                    }

                    for (int it_channel = 0; it_channel < 3; ++it_channel) {

#if STATISTICS_TIME
                        getBlock.tic();
#endif

                        lf.getBlock(orig4D, it_pos, dimBlock, stride_block, encoderParameters.dim_block, stride_lf,
                                    it_channel);

#if STATISTICS_TIME
                        getBlock.toc();
#endif

#if TRANSF_QUANT

#if STATISTICS_TIME
                        t.tic();
#endif
                        /*predictor.predict(orig4D, encoderParameters.dim_block, pf4D);
                        transform.dct_4d(pf4D, tf4D, dimBlock, encoderParameters.dim_block);*/

                        transform.dct_4d(orig4D, tf4D, dimBlock, encoderParameters.dim_block);

#if STATISTICS_TIME
                        t.toc();
#endif

#if TRACE_TRANSF
                        file_traceTransf <<
                                         it_channel << sep <<

                                         it_pos.x << sep <<
                                         it_pos.y << sep <<
                                         it_pos.u << sep <<
                                         it_pos.v << sep <<

                                         dimBlock.x << sep <<
                                         dimBlock.y << sep <<
                                         dimBlock.u << sep <<
                                         dimBlock.v << sep;

                        for (auto it : tf4D) {
                            file_traceTransf << it << sep;
                        }

                        file_traceTransf << std::endl;
#endif

#if QUANTIZATION

#if STATISTICS_TIME
                        q.tic();
#endif

                        quantization.foward(tf4D, qf4D);

#if DPCM_DC
                        qf4D[0] -= (float) dpcmDc[it_channel].get_reference(it_pos.x, it_pos.y);
#endif

#if STATISTICS_TIME
                        q.toc();
#endif

                        for (int i = 0; i < encoderParameters.dim_block.getNSamples(); ++i) {
                            temp_lre[i] = (int) std::trunc(qf4D[i]);
                        }

#if HEXADECA_TREE
                        root = tree.CreateRoot(light_field_name, hypercubo, it_channel, temp_lre, encoderParameters.dim_block);
                        tree.CreateTree(root, 0, it_pos, hypercubo_pos, {0,0,0,0});

#if HEXADECA_TREE_CODEC_MODE == 0
                        tree.ComputeHierarchicalCBF();
#elif HEXADECA_TREE_CODEC_MODE == 1
                        tree.ComputeLastCBF();
#elif HEXADECA_TREE_CODEC_MODE == 2
                        tree.ComputeLastRun();
#elif HEXADECA_TREE_CODEC_MODE == 3
                        tree.ComputeMultipleLevesOfLastRun();
#elif HEXADECA_TREE_CODEC_MODE == 4
                        tree.ComputeLastRun_0_1();
#endif

                        tree.DeleteTree(&root);

                        /*if (hypercubo == 1)
                           exit(1);*/
#endif
                        auto lre_result = lre.encodeCZI(temp_lre, 0, encoderParameters.dim_block.getNSamples());

                        bits_per_4D_Block = encoder.write4DBlock(temp_lre, encoderParameters.dim_block.getNSamples(), lre_result);

#if TRACE_QUANT
                        file_traceQuant <<
                                        it_channel << sep <<
                                        it_pos.x << sep <<
                                        it_pos.y << sep <<
                                        it_pos.u << sep <<
                                        it_pos.v << sep <<

                                        dimBlock.x << sep <<
                                        dimBlock.y << sep <<
                                        dimBlock.u << sep <<
                                        dimBlock.v << sep;

                        for (auto it : qf4D) {
                            file_traceQuant << (int) it << sep;
                        }

                        file_traceQuant << std::endl;
#endif

#if TRACE_LRE
                        for (auto it_lre : lre_result) {
                            file_traceLRE <<
                                          it_channel << sep <<
                                          it_pos.x << sep <<
                                          it_pos.y << sep <<
                                          it_pos.u << sep <<
                                          it_pos.v << sep <<

                                          dimBlock.x << sep <<
                                          dimBlock.y << sep <<
                                          dimBlock.u << sep <<
                                          dimBlock.v << sep <<

                                          it_lre.level << sep <<
                                          it_lre.run <<
                                          std::endl;
                        }
#endif

                        // TODO: entropy.encode

#if STATISTICS_TIME
                        qi.tic();
#endif

#if DPCM_DC
                        qf4D[0] += (float) dpcmDc[it_channel].get_reference(it_pos.x, it_pos.y);
#endif
                        quantization.inverse(qf4D, qi4D);

#if STATISTICS_TIME
                        qi.toc();
#endif

#else /* NO_QUANTIZATION */
                        std::copy(tf4D, tf4D + encoderParameters.dim_block.getNSamples(), qi4D);
#endif

#if STATISTICS_TIME
                        ti.tic();
#endif

                        transform.idct_4d(qi4D, ti4D, dimBlock, encoderParameters.dim_block);
                        //predictor.rec(ti4D, pi4D, dimBlock);

#if STATISTICS_TIME
                        ti.toc();
#endif

#else /* NO_TRANSF_QUANT */
                        std::copy(orig4D, orig4D + dimBlock.getNSamples(), ti4D);
#endif

#if STATISTICS_TIME
                        rebuild.tic();
#endif

                        //lf.rebuild(ti4D, it_pos, dimBlock, stride_block, encoderParameters.dim_block, stride_lf, it_channel);
                        lf.rebuild(pi4D, it_pos, dimBlock, stride_block, encoderParameters.dim_block, stride_lf, it_channel);

#if STATISTICS_TIME
                        rebuild.toc();
#endif

#if STATISTICS_LOCAL
#if TRANSF_QUANT
                        //                        statistics_tf.compute(std::vector<float>(tf4D, tf4D + dimBlock.getNSamples()));
                        //                        statistics_tf.write(it_pos, it_channel);
#if QUANTIZATION
                        statistics_qf.compute(
                                std::vector<float>(qf4D, qf4D + encoderParameters.dim_block.getNSamples()),
                                lre.getScanOrder());
                        statistics_qf.compute_sse(orig4D, qi4D, dimBlock, stride_block);
                        statistics_qf.write(it_pos, dimBlock, it_channel, lre_result, bits_per_4D_Block);
#endif // QUANTIZATION
#endif // TRANSF_QUANT
#endif // STATISTICS_LOCAL

#if DPCM_DC
                        dpcmDc[it_channel].update((int) qf4D[0], true);
#endif
                        encoder.write_completedBytes();
                    }
#if HEXADECA_TREE
                    hypercubo++;
#endif
                }
            }
        }
    }

    //lf.write(encoderParameters.getPathOutput());
    //encoder.finish_and_write();
    encoder.~EncBitstreamWriter();

#if STATISTICS_TIME
    total_time.toc();
#endif

    cout << "\n#########################################################" << endl;
    cout << "Total bytes:\t" << encoder.getTotalBytes() << endl;
    cout << "#########################################################" << endl;

#if HEXADECA_TREE
    tree.CloseFile();
#endif

#if STATISTICS_GLOBAL
    // TODO: statistics (global)
#endif

#if TRACE_QUANT
    file_traceQuant.close();
#endif

#if TRACE_LRE
    file_traceLRE.close();
#endif

#if STATISTICS_TIME
    std::ofstream file_time;
    file_time.open(encoderParameters.getPathOutput() + "time.csv");
    file_time <<
              "GetBlock Time(sec)" << sep << "T Time(sec)" << sep << "Q Time(sec)" << sep << "TI Time(sec)" << sep <<
              "QI Time(sec)" << sep << "Rebuild Time(sec)" << sep << "Total Time(sec)" << sep;
    file_time << std::endl;
    file_time <<
              getBlock.getTotalTime() << sep << t.getTotalTime() << sep << q.getTotalTime() << sep << ti.getTotalTime()
              << sep <<
              qi.getTotalTime() << sep << rebuild.getTotalTime() << sep << total_time.getTotalTime() << sep;
    file_time << std::endl;
#endif

    return 0;
}