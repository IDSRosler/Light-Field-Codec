
#ifndef PREDICTION_H
#define PREDICTION_H

#include "Typedef.h"
#include "clip.h"
#include "Point4D.h"
#include <string>
#include <cstring>
//EDUARDO BEGIN
#include <vector>


class ValueBlockPred {
public:
    bool available{false};

    std::vector<float> block4D;

    ValueBlockPred() = default;

    ValueBlockPred(float *block4D, bool available, uint blockSize);

};
//EDUARDO END

class Prediction {

public:
    explicit Prediction();

    ~Prediction();

    LFSample *rgb[3];
    LFSample predictors[3];
    LFSample predictor;
    int l;
    int mode_Selected[1218],  sse_Selected[1218] = {0};

    uint num_elementos;

    //EDUARDO BEGIN
    uint resol_x; //tamanho

    std::vector<ValueBlockPred> pred_references;

    void init_references();

    void writeHeatMap(const std::string  output_path);

    float sseHorizontalFullBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize);
    float sseVerticalFullBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize);
  
    Prediction(uint resol_x);

    int get_reference(uint x, uint y);

    void update(float *curr, bool available, uint blockSize);

    void get_referenceL(uint x, uint y, float *out, const Point4D &origSize, bool &available);

    void get_referenceA(uint x, uint y, float *out, const Point4D &origSize, bool &available);

    void get_referenceAL(uint x, uint y, float *out, const Point4D &origSize, bool &available);

    void get_referenceAR(uint x, uint y, float *out, const Point4D &origSize, bool &available);

    void predictRef(const float *orig_input, const float *ref, const Point4D &origSize, float *out );

    void DC(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out );

    void IBC(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out );

    float sseBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize);

    float sseHorizontal(const float *orig_input, const float *prediction_input, const Point4D &origSize);

    float sseVertical(const float *orig_input, const float *prediction_input, const Point4D &origSize);

    void generateReferenceVectorHorizontal(const float *blockRef1, bool availableRef1, const float *blockRef2, bool availableRef2, const Point4D &origSize, float *out );

    void generateReferenceVectorVertical(const float *blockRef1, bool availableRef1, const float *blockRef2, bool availableRef2, const Point4D &origSize, float *out );

    void angularPredictionVector(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out, int block, float *ref);

    void angularPrediction(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out, int block, float *ref );

    float roundTowardsZero( const float value );

    void residuePred(const float *orig_input, const float *pred, const Point4D &origSize, float *out );

    void recResiduePred(const float *orig_input, const float *pred, const Point4D &origSize, float *out );

    void YCbCR2RGB(float **yCbCr, const Point4D &origSize, float **rgb, int mPGMScale);

    void YCbCR2RGBVector(float **yCbCr, const Point4D &origSize, float **rgb, int mPGMScale);

    void write(float **rgb, const Point4D &origSize, int mPGMScale, int start_t, int start_s, const std::string fileName);

    void writeVector(float **rgb, const Point4D &origSize, int mPGMScale, int start_t, int start_s, const std::string fileName);

    void WritePixelToFile(int pixelPositionInCache, float **rgb, int mPGMScale, int mNumberOfFileBytesPerPixelComponent, FILE *mViewFilePointer);

    unsigned short change_endianness_16b(unsigned short val);

    void blockGenerator(const Point4D &origSize, int mode, int mPGMScale, int start_t, int start_s, const std::string fileName);

    void recRef(const float *input, const Point4D &origSize, float *out );

    //EDUARDO END

    int mNumberOfHorizontalViews, mNumberOfVerticalViews;

    int mColumns, mLines;

    int start_t{0}, start_s{0};

    void predict(const float *orig_input, const Point4D &origSize, float *out );

    void getBlock(float *block, const Point4D &pos, const Point4D &dim_block, const Point4D &stride_block,
                  const Point4D &origSize, const Point4D &stride_lf, int channel);

    void rec(float *input, float *out, Point4D &dim_block);

};

#endif //PREDICTION_H
