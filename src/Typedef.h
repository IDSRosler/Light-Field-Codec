
#ifndef TYPEDEF_H
#define TYPEDEF_H

#define USE_YCbCr 1 // 0: no, 1: mule, 2: other

#define LFSample short

#define TRANSF_QUANT true
#define QUANTIZATION true

#define STATISTICS_LOCAL false
#define STATISTICS_GLOBAL false
#define STATISTICS_TIME false

#define TRACE_TRANSF false
#define TRACE_QUANT false
#define TRACE_LRE false

#define DPCM_DC true

#define HEXADECA_TREE true
#define HEXADECA_TREE_CODEC_MODE 4 // 0: CBF Hierarquico || 1: LAST + CBF || 2: LAST + RUN || 3: Multiple levels of LAST + RUN || 4: LAST + RUN 0 1 || 5: Other
#define HEXADECA_TREE_TYPE 0 // 0: (0s,1s,2s,n>2) | 1: heat map
#define HEXADECA_TREE_PARTITION 2 // 0: Original | 1: Order 8 | 2: Order 4
#define HYPERCUBE_POSITION 10 //Position of hypercube for heat map analyze (15x15x13x13 range 0 - 1217)
#define HYPERCUBE_CHANNEL 0 //Channel for heat map analyze

#endif //TYPEDEF_H
