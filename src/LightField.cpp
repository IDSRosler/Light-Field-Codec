
#include "LightField.h"

LightField::LightField(Point4D &dim_lf, const std::string &path, bool isLytro) {

    this->mNumberOfHorizontalViews = dim_lf.u;
    this->mNumberOfVerticalViews = dim_lf.v;
    this->numberOfViewCacheLines = dim_lf.v;
    this->fullTag[0] = '\0';

    //IDM tirado do IDM LF stats no encoder  = par.transformLength_v=15;
    //numberOfViewCacheLines = 15;

    this->start_t = isLytro ? 1 : 0;
    this->start_s = isLytro ? 1 : 0;

    this->offset = Point4D(1, dim_lf.x,
                           dim_lf.x * dim_lf.y,
                           dim_lf.x * dim_lf.y * dim_lf.u);

    for (int i = 0; i < 3; i++) {
        this->rgb[i] = new LFSample[dim_lf.getNSamples()];
        this->yCbCr[i] = new float[dim_lf.getNSamples()];
    }

    this->read(path);
}

void LightField::read(const std::string &path) {

    int cont = 0;

    for (int index_t = this->start_t; index_t < start_t + mNumberOfVerticalViews; index_t++) {
        for (int index_s = this->start_s; index_s < start_s + mNumberOfHorizontalViews; index_s++) {
            char tag[256];
            fullTag[0] = '\0';

            std::string indice_t = std::to_string(index_t);
            std::string indice_s = std::to_string(index_s);
            indice_t = std::string(3 - indice_t.length(), '0') + indice_t;
            indice_s = std::string(3 - indice_s.length(), '0') + indice_s;

            std::string name_ppm = path + indice_t + "_" + indice_s + ".ppm";

//            printf("Opening View %s \n", name_ppm.c_str());

            mViewFilePointer = fopen(name_ppm.c_str(), "rb");

            if (mViewFilePointer == nullptr) {
                printf("unable to open %s view file for reading\n", name_ppm.c_str());
                assert(false);
            }

            fscanf(mViewFilePointer, "%s", tag);

            if (strcmp(tag, "P6") == 0) {
                strcat(fullTag, tag);
                strcat(fullTag, "\n");
                fscanf(mViewFilePointer, "%s", tag);
                strcat(fullTag, tag);
                strcat(fullTag, " ");

                while (tag[0] == '#') {
                    while (tag[0] != '\n') {
                        tag[0] = fgetc(mViewFilePointer);
                        strcat(fullTag, tag);
                    }
                    fscanf(mViewFilePointer, "%s", tag);
                    strcat(fullTag, tag);
                    strcat(fullTag, "\n");
                }
                mColumns = atoi(tag);

                fscanf(mViewFilePointer, "%d", &mLines);
                strcat(fullTag, (std::to_string(mLines).c_str()));
                strcat(fullTag, "\n");

                fscanf(mViewFilePointer, "%d", &mPGMScale);
                strcat(fullTag, (std::to_string(mPGMScale).c_str()));


                fgetc(mViewFilePointer);    //reads newline character
                strcat(fullTag, "\0");

                long mPGMDataOffset = ftell(mViewFilePointer);
            } else {
                printf("Error: input file is not PPM\n");
                return;
            }

            mFirstPixelPosition = cont * mColumns * mLines;
            mNumberOfFileBytesPerPixelComponent = (mPGMScale <= 255 ? 1 : 2);
            for (int pixelCount = 0; pixelCount < mColumns * mLines; pixelCount++) {
                ReadPixelFromFile(pixelCount);
            }

            cont++;

        }
    }

    this->RGB2YCbCr();
}

void LightField::ReadPixelFromFile(int pixelPosition) {
    LFSample pixelValue;

    for (int component_index = 0; component_index < 3; component_index++) {
        if (mViewFilePointer == nullptr) {
            perror("Failed: ");
            assert(false);
        }

        fread(&pixelValue, mNumberOfFileBytesPerPixelComponent, 1, mViewFilePointer);

        this->rgb[component_index][mFirstPixelPosition + pixelPosition] =
                (mNumberOfFileBytesPerPixelComponent == 2) ? change_endianness_16b(pixelValue)
                                                           : pixelValue;
    }
}

void LightField::write(const std::string &path) {
    this->YCbCR2RGB();

    int cont = 0;
    int cont2 = 0;
    for (int index_t = this->start_t; index_t < start_t + mNumberOfVerticalViews; index_t++) {
        for (int index_s = this->start_s; index_s < start_s + mNumberOfHorizontalViews; index_s++) {

            std::string indice_t = std::to_string(index_t);
            std::string indice_s = std::to_string(index_s);
            indice_t = std::string(3 - indice_t.length(), '0') + indice_t;
            indice_s = std::string(3 - indice_s.length(), '0') + indice_s;

            std::string name_ppm = path + indice_t + "_" + indice_s + ".ppm";


            //printf("Opening View %s \n", name_ppm.c_str());
            mViewFilePointer = fopen(name_ppm.c_str(), "w");
            if (mViewFilePointer == nullptr) {
                printf("unable to open %s view file for writing\n", name_ppm.c_str());
                assert(false);
            }

            //fwrite(fullTag, sizeof(char), sizeof(fullTag), mViewFilePointer);


            mNumberOfFileBytesPerPixelComponent = (mPGMScale <= 255 ? 1 : 2);

            fprintf(mViewFilePointer, "P6\n%d %d\n%d\n", mColumns, mLines, mPGMScale);

            mFirstPixelPosition = cont * mColumns * mLines;
            for (int pixelCount = 0; pixelCount < mColumns * mLines; pixelCount++) {
                WritePixelToFile(pixelCount);
            }

            cont++;
        }
    }
}


void LightField::WritePixelToFile(int pixelPositionInCache) {

    for (int component_index = 0; component_index < 3; component_index++) {
        int ClippedPixelValue = this->rgb[component_index][mFirstPixelPosition + pixelPositionInCache];
        if (ClippedPixelValue > mPGMScale)
            ClippedPixelValue = mPGMScale;
        if (ClippedPixelValue < 0)
            ClippedPixelValue = 0;
        unsigned short bigEndianPixelValue = (mNumberOfFileBytesPerPixelComponent == 2) ? change_endianness_16b(
                ClippedPixelValue) : ClippedPixelValue;

        fwrite(&bigEndianPixelValue, mNumberOfFileBytesPerPixelComponent, 1, mViewFilePointer);
    }

}

LightField::~LightField() {

    for (int i = 0; i < 3; ++i) {
        delete[] this->rgb[i];

        delete[] this->yCbCr[i];
    }
}

unsigned short LightField::change_endianness_16b(unsigned short val) {
    return (val << 8u) | ((val >> 8u) & 0x00ff);
}


void LightField::getBlock(float *block, const Point4D &pos, const Point4D &dim_block, const Point4D &stride_block,
                          const Point4D &origSize, const Point4D &stride_lf, int channel) {
    float *it_block = block;

    float *it_input = &this->yCbCr[channel][pos.x * this->offset.x + pos.y * this->offset.y +
                                            pos.u * this->offset.u + pos.v * this->offset.v];

    for (int it_v = 0; it_v < dim_block.v; ++it_v) {
        for (int it_u = 0; it_u < dim_block.u; ++it_u) {
            for (int it_y = 0; it_y < dim_block.y; ++it_y) {
                for (int it_x = 0; it_x < dim_block.x; ++it_x) {
                    *it_block = *it_input;

                    it_block += stride_block.x;
                    it_input += stride_lf.x;
                }
                it_block += stride_block.y;
                it_input += stride_lf.y;
            }
            it_block += stride_block.u;
            it_input += stride_lf.u;
        }
        it_block += stride_block.v;
        it_input += stride_lf.v;
    }
}

void LightField::rebuild(float *block, const Point4D &pos, const Point4D &dim_block, const Point4D &stride_block,
                         const Point4D &origSize,
                         const Point4D &stride_lf, int channel) {
    float *it_block = block;

    auto *it_output = &this->yCbCr[channel][pos.x * this->offset.x + pos.y * this->offset.y +
                                            pos.u * this->offset.u + pos.v * this->offset.v];

    for (int it_v = 0; it_v < dim_block.v; ++it_v) {
        for (int it_u = 0; it_u < dim_block.u; ++it_u) {
            for (int it_y = 0; it_y < dim_block.y; ++it_y) {
                for (int it_x = 0; it_x < dim_block.x; ++it_x) {
                    *it_output = std::round(*it_block);

                    it_block += stride_block.x;
                    it_output += stride_lf.x;
                }
                it_block += stride_block.y;
                it_output += stride_lf.y;
            }
            it_block += stride_block.u;
            it_output += stride_lf.u;
        }
        it_block += stride_block.v;
        it_output += stride_lf.v;
    }
}

void LightField::RGB2YCbCr() {
    int cont = 0;

    int N = 10;
    double nd = (double) (1 << (N - 8));  // pow(2, (double)N - 8);

    double clipval = (double) (1 << N) - 1;
    //printf("\nClipVal = %f", nd);

    double M[] = {0.212600000000000, -0.114572000000000, 0.500000000000000,
                  0.715200000000000, -0.385428000000000, -0.454153000000000,
                  0.072200000000000, 0.500000000000000, -0.045847000000000,};

    double r = 0.0, g = 0.0, b = 0.0;


    for (int index_t = 0; index_t < mNumberOfVerticalViews; index_t++) {
        for (int index_s = 0; index_s < mNumberOfHorizontalViews; index_s++) {
            mFirstPixelPosition = cont * mColumns * mLines;
            cont++;

            for (int pixelCount = 0; pixelCount < mColumns * mLines; pixelCount++) {

#if USE_YCbCr == 1 /*Mule*/

                r = (double) rgb[0][mFirstPixelPosition + pixelCount] / clipval;
                g = (double) rgb[1][mFirstPixelPosition + pixelCount] / clipval;
                b = (double) rgb[2][mFirstPixelPosition + pixelCount] / clipval;

                for (int icomp = 0; icomp < 3; icomp++) {

                    yCbCr[icomp][mFirstPixelPosition + pixelCount] = r * M[icomp + 0]
                                                                     + g * M[icomp + 3]
                                                                     + b * M[icomp + 6];

                    if (icomp < 1) {
                        yCbCr[icomp][mFirstPixelPosition + pixelCount] = (219
                                                                          * yCbCr[icomp][mFirstPixelPosition +
                                                                                         pixelCount] + 16) * nd;
                    } else {
                        yCbCr[icomp][mFirstPixelPosition + pixelCount] = (224
                                                                          * yCbCr[icomp][mFirstPixelPosition +
                                                                                         pixelCount] + 128) * nd;
                    }
                    yCbCr[icomp][mFirstPixelPosition + pixelCount] =
                            yCbCr[icomp][mFirstPixelPosition + pixelCount] - (mPGMScale + 1) / 2;
                }


#elif USE_YCbCr == 2 /*other*/

                yCbCr[0][mFirstPixelPosition + pixelCount] =
                        (((float) rgb[0][mFirstPixelPosition + pixelCount] * 65.738f / 256)) +
                        (((float) rgb[1][mFirstPixelPosition + pixelCount] * 129.057f / 256)) +
                        (((float) rgb[2][mFirstPixelPosition + pixelCount] * 25.064f / 256)) + 16.0f;

                yCbCr[1][mFirstPixelPosition + pixelCount] =
                        (-((float) rgb[0][mFirstPixelPosition + pixelCount] * 37.945f / 256)) -
                        (((float) rgb[1][mFirstPixelPosition + pixelCount] * 74.494f / 256)) +
                        (((float) rgb[2][mFirstPixelPosition + pixelCount] * 112.439f / 256)) + 128.0f;

                yCbCr[2][mFirstPixelPosition + pixelCount] =
                        (((float) rgb[0][mFirstPixelPosition + pixelCount] * 112.439f / 256)) -
                        (((float) rgb[1][mFirstPixelPosition + pixelCount] * 94.154f / 256)) -
                        (((float) rgb[2][mFirstPixelPosition + pixelCount] * 18.285f / 256)) + 128.0f;

                yCbCr[0][mFirstPixelPosition + pixelCount] =
                        yCbCr[0][mFirstPixelPosition + pixelCount] - (mPGMScale + 1) / 2;
                yCbCr[1][mFirstPixelPosition + pixelCount] =
                        yCbCr[1][mFirstPixelPosition + pixelCount] - (mPGMScale + 1) / 2;
                yCbCr[2][mFirstPixelPosition + pixelCount] =
                        yCbCr[2][mFirstPixelPosition + pixelCount] - (mPGMScale + 1) / 2;


#else /*No - RGB*/
                yCbCr[0][mFirstPixelPosition + pixelCount] =
                        (LFSample) rgb[0][mFirstPixelPosition + pixelCount];
                yCbCr[1][mFirstPixelPosition + pixelCount] =
                        (LFSample) rgb[1][mFirstPixelPosition + pixelCount];
                yCbCr[2][mFirstPixelPosition + pixelCount] =
                        (LFSample) rgb[2][mFirstPixelPosition + pixelCount];
#endif
            }

        }
    }

}

void LightField::YCbCR2RGB() {

    int cont = 0;

    int N = 10;
    float pixel[3];
    double M[] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 0,
                  -0.187330000000000, 1.855630000000000, 1.574800000000000,
                  -0.468130000000000, 0};


    double nd = (double) (1 << (N - 8));

    unsigned short clipval = (unsigned short) (1 << N) - 1;  // pow(2, N) - 1;

    double sval1 = 16 * nd;
    double sval2 = 219 * nd;
    double sval3 = 128 * nd;
    double sval4 = 224 * nd;


    for (int index_t = 0; index_t < mNumberOfVerticalViews; index_t++) {
        for (int index_s = 0; index_s < mNumberOfHorizontalViews; index_s++) {
            mFirstPixelPosition = cont * mColumns * mLines;
            cont++;

            for (int pixelCount = 0; pixelCount < mColumns * mLines; pixelCount++) {

#if USE_YCbCr == 1 /*Mule*/

                for (int icomp = 0; icomp < 3; icomp++) {
                    yCbCr[icomp][mFirstPixelPosition + pixelCount] =
                            yCbCr[icomp][mFirstPixelPosition + pixelCount] + (mPGMScale + 1) / 2;

                    if (icomp < 1) {
                        yCbCr[icomp][mFirstPixelPosition + pixelCount] = clip(
                                (yCbCr[icomp][mFirstPixelPosition + pixelCount] - sval1) / sval2, 0.0, 1.0);
                    } else {
                        yCbCr[icomp][mFirstPixelPosition + pixelCount] = clip(
                                (yCbCr[icomp][mFirstPixelPosition + pixelCount] - sval3) / sval4, -0.5, 0.5);
                    }

                }

                for (int icomp = 0; icomp < 3; icomp++) {

                    pixel[icomp] = yCbCr[0][mFirstPixelPosition + pixelCount] * M[icomp + 0]
                                   + yCbCr[1][mFirstPixelPosition + pixelCount] * M[icomp + 3]
                                   + yCbCr[2][mFirstPixelPosition + pixelCount] * M[icomp + 6];

                    rgb[icomp][mFirstPixelPosition + pixelCount] = std::round(clip(
                            double(pixel[icomp] * clipval), 0.0, (double) clipval));
                }


#elif USE_YCbCr == 2 /*Other*/
                yCbCr[0][mFirstPixelPosition+pixelCount] = yCbCr[0][mFirstPixelPosition+pixelCount] + (mPGMScale+1)/2;
                yCbCr[1][mFirstPixelPosition+pixelCount] = yCbCr[1][mFirstPixelPosition+pixelCount] + (mPGMScale+1)/2;
                yCbCr[2][mFirstPixelPosition+pixelCount] = yCbCr[2][mFirstPixelPosition+pixelCount] + (mPGMScale+1)/2;

                rgb[0][mFirstPixelPosition + pixelCount] = (LFSample) round(
                        (298.082f / 256 * (yCbCr[0][mFirstPixelPosition + pixelCount] - 16.0f)) +
                        (408.583f / 256 * (yCbCr[2][mFirstPixelPosition + pixelCount] - 128.0f)));


                rgb[1][mFirstPixelPosition + pixelCount] = (LFSample) round(
                        (298.082f / 256 * (yCbCr[0][mFirstPixelPosition + pixelCount] - 16.0f)) -
                        (208.120f / 256 * (yCbCr[2][mFirstPixelPosition + pixelCount] - 128.0f)) -
                        (100.291f / 256 * (yCbCr[1][mFirstPixelPosition + pixelCount] - 128.0f)));

                rgb[2][mFirstPixelPosition + pixelCount] = (LFSample) round(
                        (298.082f / 256 * (yCbCr[0][mFirstPixelPosition + pixelCount] - 16.0f)) +
                        (516.412f / 256 * (yCbCr[1][mFirstPixelPosition + pixelCount] - 128.0f)));


#else /*No - RGB*/
                rgb[0][mFirstPixelPosition + pixelCount] =
                        (LFSample) yCbCr[0][mFirstPixelPosition + pixelCount];
                rgb[1][mFirstPixelPosition + pixelCount] =
                        (LFSample) yCbCr[1][mFirstPixelPosition + pixelCount];
                rgb[2][mFirstPixelPosition + pixelCount] =
                        (LFSample) yCbCr[2][mFirstPixelPosition + pixelCount];
#endif
            }
        }

    }

}
