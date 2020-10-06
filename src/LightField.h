
#ifndef LIGHTFIELD_H
#define LIGHTFIELD_H

#include "Typedef.h"
#include "clip.h"
#include "Point4D.h"
#include <string>
#include <cmath>
#include <cstring>
#include <cassert>


class LightField {
    void read(const std::string &path);

    inline unsigned short change_endianness_16b(unsigned short val);

public:
    explicit LightField(Point4D &dim_lf, const std::string &path, bool isLytro);

    ~LightField();

    LFSample *rgb[3];
    float *yCbCr[3];

    Point4D offset;

    void write(const std::string &path);

    int mNumberOfHorizontalViews, mNumberOfVerticalViews;
    int start_t{0}, start_s{0};
    int mFirstPixelPosition, mNumberOfFileBytesPerPixelComponent;

    char fullTag[265];
    FILE *mViewFilePointer;
    int mColumns, mLines;
    int mPGMScale;

    void ReadPixelFromFile(int pixelPosition);

    void WritePixelToFile(int pixelPositionInCache);

    int numberOfViewCacheLines;

    void
    getBlock(float *block, const Point4D &pos, const Point4D &dim_block, const Point4D &stride_block,
             const Point4D &origSize, const Point4D &stride_lf, int channel);

    void
    rebuild(float *block, const Point4D &pos, const Point4D &dim_block, const Point4D &stride_block,
            const Point4D &origSize, const Point4D &stride_lf, int channel);

    void RGB2YCbCr();

    void YCbCR2RGB();
};


#endif //LIGHTFIELD_H
