
#ifndef TYPEDEF_H
#define TYPEDEF_H

#define USE_YCbCr 1 // 0: no, 1: mule, 2: other

using index_t = unsigned;

#define LFSample short



#ifndef DPCM_DC
#define DPCM_DC false
#endif
#ifndef STATISTICS_LOCAL
#define STATISTICS_LOCAL false
#endif
#ifndef STATISTICS_GLOBAL
#define STATISTICS_GLOBAL false
#endif
#ifndef STATISTICS_TIME
#define STATISTICS_TIME true
#endif
#ifndef TRACE_TRANSF
#define TRACE_TRANSF false
#endif
#ifndef TRACE_QUANT
#define TRACE_QUANT false
#endif
#ifndef TRACE_LRE
#define TRACE_LRE false
#endif
#ifndef LFCODEC_USE_PREDICTION
#define LFCODEC_USE_PREDICTION true
#endif
#ifndef LFCODEC_ESTIMATE_SCAN_ORDER
#define LFCODEC_ESTIMATE_SCAN_ORDER true
#endif

#ifndef LFCODEC_FORCE_DCT_NON_LUMA
#define LFCODEC_FORCE_DCT_NON_LUMA true
#endif

#ifndef LFCODEC_TRACE_TRANSFORM
#define LFCODEC_TRACE_TRANSFORM false
#endif

#define ENTROPY_TYPE true  // true: Arithmetic | false: LRE
#define HEXADECA_TREE_PARTITION 2 // 0: Original | 1: Order 8 | 2: Order 4

#ifndef LFCODEC_USE_QUANTIZATION
#define LFCODEC_USE_QUANTIZATION true
#endif



#define LFCODEC_EXPORT_MICROIMAGES_TRANSFORM false




#endif // TYPEDEF_H
