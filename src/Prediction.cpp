#include "Prediction.h"
#include <cmath>
#include <iostream>

//idm bibliotecas pra escrita do arquivo CSV
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

Prediction::Prediction(){
}

ValueBlockPred::ValueBlockPred(float *block4D, bool available, uint blockSize){
    this->available = available;
    for(int i = 0; i < blockSize ; ++i){
        this->block4D.push_back(block4D[i]);
    }
}

Prediction::Prediction(uint resol_x) : resol_x(resol_x) {
    this->init_references();
}

void Prediction::get_referenceL(uint x, uint y, float *out, const Point4D &origSize, bool &available) {
    ValueBlockPred ref = *(this->pred_references.end() - 1); // left
    int numElements = origSize.getNSamples();
    if (x == 0) ref.available = false;
    if(ref.available){
        for (int i = 0; i < numElements ; ++i)
            out[i] = ref.block4D[i];
        available = true;
    }else{
        for (int i = 0; i < numElements ; ++i)
            out[i] = 0;
        available = false;
    }
}


void Prediction::get_referenceAL(uint x, uint y, float *out, const Point4D &origSize, bool &available) {
    ValueBlockPred ref = *this->pred_references.begin(); // left above
    int numElements = origSize.getNSamples();
    if (y == 0) ref.available = false;
    if (x == 0) ref.available = false;
    if(ref.available){
        for (int i = 0; i < numElements ; ++i)
            out[i] = ref.block4D[i];
        available = true;
    }else{
        for (int i = 0; i < numElements ; ++i)
            out[i] = 0;
        available = false;
    }
}


void Prediction::get_referenceA(uint x, uint y, float *out, const Point4D &origSize, bool &available) {
    ValueBlockPred ref = *(this->pred_references.begin() + 1); // above
    int numElements = origSize.getNSamples();

    if (y == 0 ) {
        ref.available = false;
    }else{
        ref.available = true;
    }

    if(ref.available){
        for (int i = 0; i < numElements ; ++i)
            out[i] = ref.block4D[i];
        available = true;
    }else{
        for (int i = 0; i < numElements ; ++i)
            out[i] = 0;
        available = false;
    }
}
void Prediction::get_referenceAR(uint x, uint y, float *out, const Point4D &origSize, bool &available, int block) {
    ValueBlockPred ref = *(this->pred_references.begin() + 2); // above right
    int numElements = origSize.getNSamples();
    this->countLFend;

    if (y == 0 ) {
        ref.available = false;
    }else{
        ref.available = true;
    }

    if (block == countLFend) {
        ref.available = false;
        countLFend += 42;
        // std::cout << countLFend << std::endl;
    }
    

    if(ref.available){
        for (int i = 0; i < numElements ; ++i)
            out[i] = ref.block4D[i];
        available = true;
    }else{
        for (int i = 0; i < numElements ; ++i)
            out[i] = 0;
        available = false;
    }
}

void Prediction::DC(uint pos_x, uint pos_y, int block, const float *orig_input, const Point4D &origSize, float *out, int channel){
    float refAbove4D[origSize.getNSamples()],
            refLeft4D[origSize.getNSamples()],
            refAboveRight4D[origSize.getNSamples()],
            refAboveLeft4D[origSize.getNSamples()];

    bool availableL, availableA, availableAR, availableAL;
    float medL = 0, medA = 0, medAR = 0,  medAL = 0, medTotal = 0;
    int cont = 0;

    this->get_referenceA(pos_x, pos_y, refAbove4D, origSize, availableA);
    this->get_referenceL(pos_x, pos_y, refLeft4D, origSize, availableL);
    this->get_referenceAR(pos_x, pos_y, refAboveRight4D, origSize, availableAR, block);
    this->get_referenceAL(pos_x, pos_y, refAboveLeft4D, origSize, availableAL);

    if(availableL){
        for(int i = 0; i < origSize.getNSamples(); i++){
            medL += refLeft4D[i];
        }
        medL = medL/origSize.getNSamples();
        medTotal += medL;
        cont++;
    }
    if(availableA){
        for(int i = 0; i < origSize.getNSamples(); i++){
            medA += refAbove4D[i];
        }
        medA = medA/origSize.getNSamples();
        medTotal += medA;
        cont++;
    }
    if(availableAR){
        for(int i = 0; i < origSize.getNSamples(); ++i){
            medAR += refAboveRight4D[i];
        }
        medAR = medAR/origSize.getNSamples();
        medTotal += medAR;
        cont++;
    }
    if(availableAL){
        for(int i = 0; i < origSize.getNSamples(); ++i){
            medAL += refAboveLeft4D[i];
        }
        medAL = medAL/origSize.getNSamples();
        medTotal += medAL;
        cont++;
    }
    if(cont == 0){
        medTotal = 0;
    } else{
        medTotal = medTotal/cont;
    }
    //for (int i = 0; i < origSize.getNSamples(); ++i){
    //    out[i] = medTotal;
    //}

    Point4D it_pos;


    for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {
        for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {
            for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
                for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

                    int pos_out = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                                  + (it_pos.v * origSize.x * origSize.y * origSize.u);

                    if (it_pos.u == 0 && (it_pos.v == 0 || it_pos.v == origSize.v - 1)
                        || it_pos.u == origSize.u - 1 && (it_pos.v == 0 || it_pos.v == origSize.v - 1)) {//pixels pretos

                        if(channel == 0){
                            out[pos_out] = -448;
                        }else{
                            out[pos_out] = 0;
                        }

                    } else {
                        out[pos_out] = medTotal;
                    }
                }
            }
        }
    }
}

void Prediction::IBC(uint pos_x, uint pos_y, int block, const float *orig_input, const Point4D &origSize, float *out ){
    float refAbove4D[origSize.getNSamples()],
            refLeft4D[origSize.getNSamples()],
            refAboveRight4D[origSize.getNSamples()],
            refAboveLeft4D[origSize.getNSamples()];

    bool availableL, availableA, availableAR, availableAL;
    float sseA = 0, sseL = 0, sseAR = 0, sseAL = 0, sse = 0;
    int index = 0;

    this->get_referenceA(pos_x, pos_y, refAbove4D, origSize, availableA);
    this->get_referenceL(pos_x, pos_y, refLeft4D, origSize, availableL);
    this->get_referenceAR(pos_x, pos_y, refAboveRight4D, origSize, availableAR, block);
    this->get_referenceAL(pos_x, pos_y, refAboveLeft4D, origSize, availableAL);

    if(availableA) {
        sseA = this->sseBlock(orig_input, refAbove4D, origSize);
        sse = sseA;
        index = 1;
    }
    if(availableL) {
        sseL = this->sseBlock(orig_input, refLeft4D, origSize);
        if(availableA){
            if(sseL < sse){
                sse = sseL;
                index = 2;
            }
        } else{
            sse = sseL;
        }
    }
    if(availableAR) {
        sseAR = this->sseBlock(orig_input, refAboveRight4D, origSize);
        if(availableA || availableL){
            if(sseAR < sse){
                sse = sseAR;
                index = 3;
            }
        } else{
            sse = sseAR;
        }
    }
    if(availableAL) {
        sseAL = this->sseBlock(orig_input, refAboveLeft4D, origSize);
        if(availableA || availableL || availableAR){
            if(sseAL < sse){
                sse = sseAL;
                index = 4;
            }
        } else{
            sse = sseAR;
        }
    }

    if(index == 1){
        for (int i = 0; i < origSize.getNSamples(); ++i)
            out[i] = refAbove4D[i];
    } else if(index == 2){
        for (int i = 0; i < origSize.getNSamples(); ++i)
            out[i] = refLeft4D[i];
    } else if(index == 3){
        for (int i = 0; i < origSize.getNSamples(); ++i)
            out[i] = refAboveRight4D[i];
    } else if(index == 4){
        for (int i = 0; i < origSize.getNSamples(); ++i)
            out[i] = refAboveLeft4D[i];
    } else if(index == 0){
        for (int i = 0; i < origSize.getNSamples(); ++i)
            out[i] = 0;
    }
}

float Prediction::sseBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize){
    float sum = 0;
    for (int i = 0; i < origSize.getNSamples() ; ++i)
        sum += pow(orig_input[i] - prediction_input[i], 2);
    return sum;
}

float Prediction::sseHorizontal(const float *orig_input, const float *prediction_input, const Point4D &origSize){
    Point4D it_pos;

    // Horizontal
    it_pos.x = origSize.x - 1;
    it_pos.u = floor(origSize.u / 2);

    it_pos.y = 0;
    it_pos.v = 0;

    float sum = 0;
    int pos = 0;

    for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {

        for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {

            pos = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                  + (it_pos.v * origSize.x * origSize.y * origSize.u);
            sum += pow(orig_input[pos] - prediction_input[pos], 2);
        }
    }
    return sum;
}

float Prediction::sseVertical(const float *orig_input, const float *prediction_input, const Point4D &origSize){
    Point4D it_pos;

    // Horizontal - variable
    it_pos.x = 0;
    it_pos.u = 0;

    // Vertical - fixed
    it_pos.y = origSize.y - 1;
    it_pos.v = floor(origSize.v / 2);

    float sum = 0;
    int pos = 0;

    for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {

        for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

            pos = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                  + (it_pos.v * origSize.x * origSize.y * origSize.u);
            sum += pow(orig_input[pos] - prediction_input[pos], 2);
        }
    }
    return sum;
}

void Prediction::generateReferenceVectorHorizontal(const float *blockRef1, bool availableRef1, const float *blockRef2, bool availableRef2, const Point4D &origSize, float *out ){
    Point4D it_pos_out;
    Point4D it_pos_in;
    it_pos_out.x = origSize.x - 1; // fixed
    it_pos_out.y = 0;
    it_pos_out.u = 0;
    it_pos_out.v = 0;

    int cont = 0, cont2 = origSize.y * origSize.u * origSize.v; //13*13*15 = 2535
    int pos_ref, pos_out;
    availableRef2 = false;

    for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) { // spatial
        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

            for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) { // angular

                pos_ref = it_pos_out.x + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.y * origSize.x) +
                          (it_pos_out.v * origSize.u * origSize.y * origSize.x);

                pos_out = it_pos_out.y + (it_pos_out.u * origSize.y) + (it_pos_out.v * origSize.u * origSize.y);

                out[pos_out] = blockRef1[pos_ref];

                if (availableRef2) {
                    out[cont2 + pos_out] = blockRef2[pos_ref];
                } else {
                    out[cont2 + pos_out] = blockRef1[pos_ref];
                }
            }
        }
    }
}

void Prediction::generateReferenceVectorVertical(const float *blockRef1, bool availableRef1, const float *blockRef2, bool availableRef2, const Point4D &origSize, float *out ){
    Point4D it_pos_out;
    Point4D it_pos_in;
    it_pos_out.x = 0;
    it_pos_out.y = origSize.y - 1; // fixed
    it_pos_out.u = 0;
    it_pos_out.v = 0;

    int cont = 0, cont2 = origSize.x * origSize.u * origSize.v; //15*13*13 = 2535
    int pos_ref, pos_out;

    for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) { // spatial
        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

            for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) { // angular

                pos_ref = it_pos_out.x + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.y * origSize.x) + (it_pos_out.v * origSize.u * origSize.y * origSize.x);

                pos_out = it_pos_out.x + (it_pos_out.u * origSize.x) + (it_pos_out.v * origSize.x * origSize.u);

                out[pos_out] = blockRef1[pos_ref];

                if(availableRef2){
                    out[cont2 + pos_out] = blockRef2[pos_ref];
                }else{
                    out[cont2 + pos_out] = blockRef1[pos_ref];
                }
            }
        }
    }
}


// void Prediction::fillRefVertical(const float *blockRef1, bool availableRef1, const float *blockRef2, bool availableRef2, const Point4D &origSize, float *out ){
 

// }


void Prediction::angularPredictionVector(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out, int block, float *ref, int channel, int mPGMScale, std::string outputPath){
    Point4D it_pos_in;
    Point4D it_pos_out;

    it_pos_in.x = 0;
    it_pos_in.u = 0;
    it_pos_in.v = 0;
    it_pos_in.y = 0;

    it_pos_out.x = 0;
    it_pos_out.u = 0;
    it_pos_out.v = 0;
    it_pos_out.y = 0;

    int num_modes = 0; //33
    int d = 0;
    int C = 0;
    int ind = 0;
    int pos = 0;
    int W = 0;
    float R0 = 0;
    float R1 = 0;
    float min_sse = 0;
    float min_mode = 0;
    int min_d = 0;

    float refAbove4D[origSize.getNSamples()],
            refLeft4D[origSize.getNSamples()],
            refAboveRight4D[origSize.getNSamples()];

    float refAboveGeneratedVector[(origSize.x * origSize.u * origSize.v)*2];
    float refLeftGeneratedVector[(origSize.y * origSize.u * origSize.v)*2];

    bool availableL, availableA, availableAR;

    this->get_referenceA(pos_x, pos_y, refAbove4D, origSize, availableA);
    this->get_referenceL(pos_x, pos_y, refLeft4D, origSize, availableL);
    this->get_referenceAR(pos_x, pos_y, refAboveRight4D, origSize, availableAR, block);

     /*std::cout << availableL << " - " << availableA << " - " << availableAR << std::endl ;*/

//IDM CHECAR CADA CASO DE IF DEPOIS
    if (not availableL && availableA){
        std::cout << "Not L and A" << std::endl;
        for(int i = 0; i < origSize.getNSamples(); i++){
                out[i] = refAbove4D[i];
        }
        //IDM REKAME PARA 512
    }else if(not availableL && not availableA){
//        std::cout << "Not L and Not A" << std::endl;
        for(int i = 0; i < origSize.getNSamples(); i++)
                out[i] = orig_input[i];
        
    }else if(not availableL ^ not availableA){ //sem uma das referências
       if(availableL){
//           std::cout << "L and Not A" << std::endl;
           for(int i = 0; i < origSize.getNSamples(); i++)
               out[i] = refLeft4D[i];
       }else{
//           std::cout << "Not L and A" << std::endl;
           for(int i = 0; i < origSize.getNSamples(); i++)
               out[i] = refAbove4D[i];
       }
            // std::cout << " ENTROU NO 1" << std::endl ;

    }else if( not availableAR ){{
//        std::cout << "Not AR" << std::endl;
        for(int i = 0; i < origSize.getNSamples(); i++)
            out[i] = refLeft4D[i];  
        
        // std::cout << " ENTROU NO 2" << std::endl ;
    }
    }else{ //com referência
        if(availableA){
//            std::cout << "Entrou em (com referência) availableA" << std::endl ;
            this->generateReferenceVectorVertical(refAbove4D, availableA, refAboveRight4D, availableAR, origSize, refAboveGeneratedVector);
        }
        if(availableL){
//            std::cout << "Entrou em (com referência) availableL" << std::endl ;
            this->generateReferenceVectorHorizontal(refLeft4D, availableL, refLeft4D, availableL, origSize, refLeftGeneratedVector);
        }

        for(int mode = 0; mode < num_modes; mode++) { //num_modes
            switch (mode)
            {
                case 0:
                case 32:
                    d = 32;
                    break;
                case 1:
                case 31:
                    d = 26;
                    break;
                case 2:
                case 30:
                    d = 21;
                    break;
                case 3:
                case 29:
                    d = 17;
                    break;
                case 4:
                case 28:
                    d = 13;
                    break;
                case 5:
                case 27:
                    d = 9;
                    break;
                case 6:
                case 26:
                    d = 5;
                    break;
                case 7:
                case 25:
                    d = 2;
                    break;
                case 8:
                case 24:
                    d = 0;
                    break;
                case 9:
                case 23:
                    d = -2;
                    break;
                case 10:
                case 22:
                    d = -5;
                    break;
                case 11:
                case 21:
                    d = -9;
                    break;
                case 12:
                case 20:
                    d = -13;
                    break;
                case 13:
                case 19:
                    d = -17;
                    break;
                case 14:
                case 18:
                    d = -21;
                    break;
                case 15:
                case 17:
                    d = -26;
                    break;
                case 16:
                    d = -32;
                    break;
                default:
                    d = 0;
            }

            if(mode <= 15 && availableL) { //Horizontal

                // Horizontal - fixed
                it_pos_in.x = origSize.x - 1;
                it_pos_in.u = floor(origSize.u / 2) * origSize.x * origSize.y;

                // Vertical - variable
                it_pos_in.v = 0;
                it_pos_in.y = 0;

                // percorre vetor out na ordem horizontal espacial
                for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {
                    for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {

                        // percorre vetor out na ordem horizontal angular
                        for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {
                            for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

                                int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) +
                                              (it_pos_out.u * origSize.x * origSize.y)
                                              + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                                //C = (it_pos_out.u * d) >> 5;
                                C = (int)this->roundTowardsZero((int)(it_pos_out.u * d) / (float)pow(2, 5));
                                W = (it_pos_out.u * d) & 31;
                                ind = it_pos_out.v + C;
                                pos = ind + 1;

                                if(ind < 0){
                                    ind = 0;
                                }

                                if (pos >= origSize.v) {
                                    pos = ind;
                                }
                                R0 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_in.u) +
                                               (ind * origSize.x * origSize.y * origSize.u)];
                                R1 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_in.u) +
                                               ((pos) * origSize.x * origSize.y * origSize.u)];

                                out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                            }
                        }
                    }
                }
                int sse = this->sseHorizontal(orig_input, out, origSize);
                if(mode == 0){
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                } else if (sse < min_sse) {
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                }

            } else if(mode > 15 && availableA){ //Vertical

                // Horizontal - variable
                it_pos_in.x = 0;
                it_pos_in.u = 0;

                // Vertical - fixed
                it_pos_in.y = (origSize.y - 1) * origSize.x;
                it_pos_in.v = floor(origSize.v / 2) * origSize.x * origSize.y * origSize.u;

                // percorre vetor out na ordem vertical espacial
                for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {
                    for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {

                        // percorre vetor out na ordem vertical angular
                        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {
                            for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {

                                int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                              + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                                //C = (it_pos_out.v * d) >> 5;
                                C = (int)this->roundTowardsZero((int)(it_pos_out.v * d) / (float)pow(2, 5));
                                W = (it_pos_out.v * d) & 31;
                                ind = it_pos_out.u + C;
                                pos = ind + 1;

                                if(ind < 0){
                                    ind = 0;
                                }

                                if (pos >= origSize.v){
                                    pos = ind;
                                }
                                R0 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + (ind * origSize.x * origSize.y) +
                                                (it_pos_in.v)];
                                R1 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + ((pos) * origSize.x * origSize.y) +
                                                (it_pos_in.v)];

                                out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                            }
                        }
                    }
                }
                int sse = this->sseVertical(orig_input, out, origSize);
                if(mode == 16 && not availableL){ //primeiro e não passou pelo horizontal
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                } else if (sse < min_sse) {
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                }

            }
        }

        //if(block == 0){
        //    std::cout << "mode: " << min_mode + 2 << " sse: " << min_sse << " d: " << min_d << std::endl;
        //}

        min_mode = 16; //fix mode
        min_d = 0; //fix d

        if(min_mode <= 15 ){ //Horizontal
            //Reaproveitamento da funcao pra preencher o vetor de referencias para escrita em arquivo

/*            std::cout << "Horizontal mode" << std::endl ;

            std::cout << "Left: " << availableL << std::endl ;
            std::cout << "Above: " << availableA << std::endl ;
            std::cout << "Above_right: " << availableAR << std::endl ;*/

            this->generateReferenceVectorHorizontal(refLeft4D, availableL, refLeft4D, availableL, origSize, ref);

            // Horizontal - fixed
            it_pos_in.x = origSize.x - 1;
            it_pos_in.u = floor(origSize.u / 2) * origSize.x * origSize.y;

            // Vertical - variable
            it_pos_in.v = 0;
            it_pos_in.y = 0;

            // percorre vetor out na ordem horizontal espacial
            for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {
                for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {

                    // percorre vetor out na ordem horizontal angular
                    for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {
                        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

                            int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                          + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                            if (it_pos_out.u == 0 && (it_pos_out.v == 0 || it_pos_out.v == origSize.v - 1)
                                || it_pos_out.u == origSize.u - 1 && (it_pos_out.v == 0 || it_pos_out.v == origSize.v - 1)) {//pixels pretos

                                if(channel == 0){
                                    out[pos_out] = -448;
                                }else{
                                    out[pos_out] = 0;
                                }

                            } else {

                                C = (int) this->roundTowardsZero((int) (it_pos_out.u * min_d) / (float) pow(2, 5));
                                W = (it_pos_out.u * min_d) & 31;
                                ind = it_pos_out.v + C;

                                if (ind < 0) {
                                    ind = 0;
                                }
                                pos = ind + 1;


//                                if(((ind == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
//                                    || (ind == 1 && it_pos_out.u <= 1)
//                                    || ((ind == 2 || ind == 3 || ind > origSize.v - 3) && it_pos_out.u == 0)) ||
//                                   ((pos == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
//                                    || (pos == 1 && it_pos_out.u <= 1)
//                                    || ((pos == 2 || pos == 3 || ind > origSize.v - 3) && it_pos_out.u == 0))){
//
//                                    out[pos_out] = refLeft4D[pos_out];
//                                } else{

                                    R0 = refLeftGeneratedVector[(it_pos_out.y) + (it_pos_out.u * origSize.y) +
                                                                (ind * origSize.y * origSize.u)];
                                    R1 = refLeftGeneratedVector[(it_pos_out.y) + (it_pos_out.u * origSize.y) +
                                                                (pos * origSize.y * origSize.u)];

                                    // std::cout << "\n Referencias: " << R0 << " - " << R1 << std::endl;

                                    /*
                                    R0 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y) +
                                                   (ind * origSize.x * origSize.y * origSize.u)];
                                    R1 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y) +
                                                   ((pos) * origSize.x * origSize.y * origSize.u)];
                                    */

                                    out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                                    
                                //}
                            }
                        }
                    }
                }
            }
        } else{ //Vertical

            // for(int i = 0; i < (origSize.x * origSize.u * origSize.v)*2; i++){
            //     ref[i] = refAboveGeneratedVector[i];
            // }

/*            std::cout << "Vertical mode" << std::endl ;

            std::cout << "Left: " << availableL << std::endl ;
            std::cout << "Above: " << availableA << std::endl ;
            std::cout << "Above_right: " << availableAR << std::endl ;*/

            this->generateReferenceVectorVertical(refAbove4D, availableA, refAboveRight4D, availableAR, origSize, ref);
            // Horizontal - variable
            it_pos_in.x = 0;
            it_pos_in.u = 0;

            // Vertical - fixed
            it_pos_in.y = (origSize.y - 1) * origSize.x;
            it_pos_in.v = floor(origSize.v / 2) * origSize.x * origSize.y * origSize.u;

            // percorre vetor out na ordem vertical espacial
            for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {
                for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {

                    // percorre vetor out na ordem vertical angular
                    for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {
                        for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {

                            int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                          + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                            if (it_pos_out.u == 0 && (it_pos_out.v == 0 || it_pos_out.v == origSize.v - 1)
                                || it_pos_out.u == origSize.u - 1 && (it_pos_out.v == 0 || it_pos_out.v == origSize.v - 1)) {//pixels pretos

                                if(channel == 0){
                                    out[pos_out] = -448;
                                }else{
                                    out[pos_out] = 0;
                                }

                            } else {

                                C = (int) this->roundTowardsZero((int) (it_pos_out.v * min_d) / (float) pow(2, 5));
                                W = (it_pos_out.v * min_d) & 31;
                                ind = it_pos_out.u + C;

                                if (ind < 0) {
                                    ind = 0;
                                }
                                pos = ind + 1;

//                                if(((ind == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
//                                    || (ind == 1 && it_pos_out.u <= 1)
//                                    || ((ind == 2 || ind == 3 || ind > origSize.v - 3) && it_pos_out.u == 0)) ||
//                                   ((pos == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
//                                    || (pos == 1 && it_pos_out.u <= 1)
//                                    || ((pos == 2 || pos == 3 || ind > origSize.v - 3) && it_pos_out.u == 0))){
//
//                                    out[pos_out] = refAbove4D[pos_out];
//                                } else {
                                    /*
                                    if(ind > origSize.u){
                                        int aux = origSize.u - ind;
                                            int contU = 1;
                                        if(aux > origSize.u){
                                            aux = origSize.u - ind;

                                        }
                                    }
                                     */

                                    // if(it_pos_out.u == 6 && it_pos_out.v == origSize.v - 1 && block == 42){
                                    //     std::cout << "ind p1: " << ind << " pos p2: " << pos << std::endl;
                                    // }

                                    R0 = refAboveGeneratedVector[(it_pos_out.x) + (ind * origSize.x) +
                                                                 (it_pos_out.v * origSize.x * origSize.u)];
                                    R1 = refAboveGeneratedVector[(it_pos_out.x) + (pos * origSize.x) +
                                                                 (it_pos_out.v * origSize.x * origSize.u)];

                                    




                                    /*
                                    R0 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + (ind * origSize.x * origSize.y) +
                                                    (it_pos_out.v * origSize.x * origSize.y * origSize.u)];
                                    R1 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + ((pos) * origSize.x * origSize.y) +
                                                    (it_pos_out.v * origSize.x * origSize.y * origSize.u)];
                                    */
                                    out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                                //}
                            }
                        }
                    }
                }
            }
        }
    }
}


void Prediction::angularPrediction(uint pos_x, uint pos_y, const float *orig_input, const Point4D &origSize, float *out, int block, float *ref){
    Point4D it_pos_in;
    Point4D it_pos_out;

    it_pos_in.x = 0;
    it_pos_in.u = 0;
    it_pos_in.v = 0;
    it_pos_in.y = 0;

    it_pos_out.x = 0;
    it_pos_out.u = 0;
    it_pos_out.v = 0;
    it_pos_out.y = 0;

    int num_modes = 33;
    int d = 0;
    int C = 0;
    int ind = 0;
    int W = 0;
    float R0 = 0;
    float R1 = 0;
    float min_sse = 0;
    float min_mode = 0;
    int min_d = 0;
    int mode = 0;

    float refAbove4D[origSize.getNSamples()],
            refLeft4D[origSize.getNSamples()];

    bool availableL, availableA, availableAR;

    this->get_referenceA(pos_x, pos_y, refAbove4D, origSize, availableA);
    this->get_referenceL(pos_x, pos_y, refLeft4D, origSize, availableL);


    int refL = 0;
    for(int i = 0; i < origSize.getNSamples(); i++){
        refL += refLeft4D[i];
    }

    int refA = 0;
    for(int i = 0; i < origSize.getNSamples(); i++){
        refA += refAbove4D[i];
    }

    if(refA == 0 && refL == 0){
        for(int i = 0; i < origSize.getNSamples(); i++){
            out[i] = orig_input[i];
        }
    } else{ //se tem bloco de referência

        for(mode = 0; mode < num_modes; mode++) { //num_modes
            switch (mode)
            {
                case 0:
                case 32:
                    d = 32;
                    break;
                case 1:
                case 31:
                    d = 26;
                    break;
                case 2:
                case 30:
                    d = 21;
                    break;
                case 3:
                case 29:
                    d = 17;
                    break;
                case 4:
                case 28:
                    d = 13;
                    break;
                case 5:
                case 27:
                    d = 9;
                    break;
                case 6:
                case 26:
                    d = 5;
                    break;
                case 7:
                case 25:
                    d = 2;
                    break;
                case 8:
                case 24:
                    d = 0;
                    break;
                case 9:
                case 23:
                    d = -2;
                    break;
                case 10:
                case 22:
                    d = -5;
                    break;
                case 11:
                case 21:
                    d = -9;
                    break;
                case 12:
                case 20:
                    d = -13;
                    break;
                case 13:
                case 19:
                    d = -17;
                    break;
                case 14:
                case 18:
                    d = -21;
                    break;
                case 15:
                case 17:
                    d = -26;
                    break;
                case 16:
                    d = -32;
                    break;
                default:
                    d = 0;
            }

            if(mode <= 15 && refL != 0) { //Horizontal

                // Horizontal - fixed
                it_pos_in.x = origSize.x - 1;
                it_pos_in.u = floor(origSize.u / 2) * origSize.x * origSize.y;

                // Vertical - variable
                it_pos_in.v = 0;
                it_pos_in.y = 0;

                // percorre vetor out na ordem horizontal espacial
                for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {
                    for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {

                        // percorre vetor out na ordem horizontal angular
                        for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {
                            for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

                                int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) +
                                              (it_pos_out.u * origSize.x * origSize.y)
                                              + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                                //C = (it_pos_out.u * d) >> 5;
                                C = (int)this->roundTowardsZero((int)(it_pos_out.u * d) / (float)pow(2, 5));
                                W = (it_pos_out.u * d) & 31;
                                ind = it_pos_out.v + C;

                                if(ind < 0){
                                    ind = 0;
                                }

                                int pos = ind + 1;
                                if (pos >= origSize.v) {
                                    pos = ind;
                                }
                                R0 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_in.u) +
                                         (ind * origSize.x * origSize.y * origSize.u)];
                                R1 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_in.u) +
                                         ((pos) * origSize.x * origSize.y * origSize.u)];

                                out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                            }
                        }
                    }
                }
                int sse = this->sseHorizontal(orig_input, out, origSize);
                if(mode == 0){
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                } else if (sse < min_sse) {
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                }

            } else if(mode > 15 && refA != 0){ //Vertical

                // Horizontal - variable
                it_pos_in.x = 0;
                it_pos_in.u = 0;

                // Vertical - fixed
                it_pos_in.y = (origSize.y - 1) * origSize.x;
                it_pos_in.v = floor(origSize.v / 2) * origSize.x * origSize.y * origSize.u;

                // percorre vetor out na ordem vertical espacial
                for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {
                    for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {

                        // percorre vetor out na ordem vertical angular
                        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {
                            for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {

                                int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                              + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                                //C = (it_pos_out.v * d) >> 5;
                                C = (int)this->roundTowardsZero((int)(it_pos_out.v * d) / (float)pow(2, 5));
                                W = (it_pos_out.v * d) & 31;
                                ind = it_pos_out.u + C;

                                if(ind < 0){
                                    ind = 0;
                                }

                                int pos = ind + 1;
                                if (pos >= origSize.v){
                                    pos = ind;
                                }
                                R0 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + (ind * origSize.x * origSize.y) +
                                         (it_pos_in.v)];
                                R1 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + ((pos) * origSize.x * origSize.y) +
                                         (it_pos_in.v)];

                                out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                            }
                        }
                    }
                }




                int sse = this->sseVertical(orig_input, out, origSize);
                if(mode == 16 && refL == 0){ //primeiro e não passou pelo horizontal
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                } else if (sse < min_sse) {
                    min_sse = sse;
                    min_mode = mode;
                    min_d = d;
                }
            }
        }

        //std::cout << "mode: " << min_mode + 2 << " sse: " << min_sse << " d: " << min_d << std::endl;

        if(min_mode <= 15 ){ //Horizontal

            for(int i = 0; i < origSize.getNSamples(); i++){
                ref[i] = refLeft4D[i];
            }

            // Horizontal - fixed
            it_pos_in.x = origSize.x - 1;
            it_pos_in.u = floor(origSize.u / 2) * origSize.x * origSize.y;

            // Vertical - variable
            it_pos_in.v = 0;
            it_pos_in.y = 0;

            // percorre vetor out na ordem horizontal espacial
            for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {
                for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {

                    // percorre vetor out na ordem horizontal angular
                    for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {
                        for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

                            int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                          + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                            if((it_pos_out.v == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                               || (it_pos_out.v == 1 && it_pos_out.u <= 1)
                               || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && it_pos_out.u == 0)){ //pixels pretos

                                out[pos_out] = refLeft4D[pos_out];

                            } else {

                                C = (int) this->roundTowardsZero((int) (it_pos_out.u * min_d) / (float) pow(2, 5));
                                W = (it_pos_out.u * min_d) & 31;
                                ind = it_pos_out.v + C;

                                if (ind < 0) {
                                    ind = 0;
                                }

                                if(ind > origSize.v-1){
                                    ind = origSize.v-1;
                                }

                                int pos = ind + 1;
                                if (pos >= origSize.v-1) {
                                    pos = ind;
                                }

                                if(((ind == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                                   || (ind == 1 && it_pos_out.u <= 1)
                                   || ((ind == 2 || ind == 3 || ind > origSize.v - 3) && it_pos_out.u == 0)) ||
                                        ((pos == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                                         || (pos == 1 && it_pos_out.u <= 1)
                                         || ((pos == 2 || pos == 3 || ind > origSize.v - 3) && it_pos_out.u == 0))){

                                    out[pos_out] = refLeft4D[pos_out];

                                } else{
                                    R0 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y) +
                                                   (ind * origSize.x * origSize.y * origSize.u)];
                                    R1 = refLeft4D[(it_pos_in.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y) +
                                                   ((pos) * origSize.x * origSize.y * origSize.u)];

                                    out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);

                                }
                            }
                        }
                    }
                }
            }
        } else{ //Vertical

            for(int i = 0; i < origSize.getNSamples(); i++){
                ref[i] = refAbove4D[i];
            }

            // Horizontal - variable
            it_pos_in.x = 0;
            it_pos_in.u = 0;

            // Vertical - fixed
            it_pos_in.y = (origSize.y - 1) * origSize.x;
            it_pos_in.v = floor(origSize.v / 2) * origSize.x * origSize.y * origSize.u;

            // percorre vetor out na ordem vertical espacial
            for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {
                for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {

                    // percorre vetor out na ordem vertical angular
                    for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {
                        for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {

                            int pos_out = (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                          + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                            if((it_pos_out.v == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                               || (it_pos_out.v == 1 && it_pos_out.u <= 1)
                               || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && it_pos_out.u == 0)){ //pixels pretos

                                out[pos_out] = refAbove4D[pos_out];

                            } else {

                                C = (int) this->roundTowardsZero((int) (it_pos_out.v * min_d) / (float) pow(2, 5));
                                W = (it_pos_out.v * min_d) & 31;
                                ind = it_pos_out.u + C;

                                if (ind < 0) {
                                    ind = 0;
                                }

                                if(ind > origSize.u-1){
                                    ind = origSize.u-1;
                                }

                                int pos = ind + 1;
                                if (pos >= origSize.u-1) {
                                    pos = ind;
                                }

                                if(((it_pos_out.v == 0 && (ind <= 3 || ind > origSize.u - 3))
                                    || (it_pos_out.v == 1 && ind <= 1)
                                    || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && ind == 0)) ||
                                   ((it_pos_out.v == 0 && (pos <= 3 || pos > origSize.u - 3))
                                    || (it_pos_out.v == 1 && pos <= 1)
                                    || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && pos == 0))){

                                    out[pos_out] = refAbove4D[pos_out];

                                } else{
                                    R0 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + (ind * origSize.x * origSize.y) +
                                                    (it_pos_out.v * origSize.x * origSize.y * origSize.u)];
                                    R1 = refAbove4D[(it_pos_out.x) + (it_pos_in.y) + ((pos) * origSize.x * origSize.y) +
                                                    (it_pos_out.v * origSize.x * origSize.y * origSize.u)];

                                    out[pos_out] = ((32 - W) * R0 + W * R1 + 16) / pow(2, 5);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
//     this->sse_Selected[l] = mode-2 <= 15 ?  sseHorizontalFullBlock(orig_input, out, origSize) : sseVerticalFullBlock(orig_input, out, origSize);


// //IDM begin Heat Map for Mode Selected

     //mode_Selected[l] = min_mode;
    // std::cout <<  "modo = " << min_mode << '\n';
//     l++;
    
//IDM end

}

//IDM begin Heat Map for Mode Selected
void Prediction::writeHeatMap(const std::string output_path){
    std::ofstream file;
    int cont = 0;

    file.open (output_path + "heat_map.csv");

    for (int i = 1; i <= 1218; i++)
    {
 
        if((i%42) == 0) file << mode_Selected[i-1] << "\n";
        else file << mode_Selected[i-1] << ",";
        
    }
    
    file.close();

}

//IDM end

float Prediction::sseHorizontalFullBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize){
    Point4D it_pos;

    // Horizontal
    it_pos.x = 0;
    it_pos.u = 0;
    it_pos.y = 0;
    it_pos.v = 0;

    float sum = 0;
    int pos = 0;

    // percorre vetor out na ordem horizontal espacial
    for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {
        for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {

            // percorre vetor out na ordem horizontal angular
            for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
                for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

                    pos = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                          + (it_pos.v * origSize.x * origSize.y * origSize.u);
                    sum += pow(orig_input[pos] - prediction_input[pos], 2);
                }
            }
        }
    }
}

float Prediction::sseVerticalFullBlock(const float *orig_input, const float *prediction_input, const Point4D &origSize){
    Point4D it_pos;
    
    it_pos.x = 0;
    it_pos.u = 0;
    it_pos.y = 0;
    it_pos.v = 0;

    float sum = 0;
    int pos = 0;

    // percorre vetor out na ordem vertical espacial
    for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {
        for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {

            // percorre vetor out na ordem vertical angular
            for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {
                for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
                    pos = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                          + (it_pos.v * origSize.x * origSize.y * origSize.u);
                    sum += pow(orig_input[pos] - prediction_input[pos], 2);
                }
            }
        }
    }
    return sum;
}



float Prediction::roundTowardsZero( const float value ){
    float result = std::floor( std::fabs( value ) );
    return (value < 0.0) ? -result : result;
}

void Prediction::residuePred(const float *orig_input, const float *pred, const Point4D &origSize, float *out ){
    for (int i = 0; i < origSize.getNSamples(); ++i){
        out[i] = orig_input[i] - pred[i];
    }
}

void Prediction::recResiduePred(const float *orig_input, const float *pred, const Point4D &origSize, float *out ){
    for (int i = 0; i < origSize.getNSamples(); ++i){
        out[i] = orig_input[i] + pred[i];
    }
}

void Prediction::YCbCR2RGB(float **yCbCr, const Point4D &origSize, float **rgb, int mPGMScale) {

    int cont = 0;
    int mFirstPixelPosition = 0;

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


    for (int index_t = 0; index_t < origSize.v; index_t++) { //vertical angular
        for (int index_s = 0; index_s < origSize.u; index_s++) { //horizontal angular
            mFirstPixelPosition = cont * origSize.x * origSize.y;
            cont++;

            for (int pixelCount = 0; pixelCount < origSize.x * origSize.y; pixelCount++) {

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

                    rgb[icomp][mFirstPixelPosition + pixelCount] = clip(
                            double(pixel[icomp] * clipval), 0.0, (double) clipval);
                }
            }
        }
    }
}
//original
void Prediction::YCbCR2RGBVector(float **yCbCr, const Point4D &origSize, float **rgb, int mPGMScale) {

    int cont = 0;
    int mFirstPixelPosition = 0;

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


    for (int index_t = 0; index_t < origSize.v; index_t++) { //vertical angular
        for (int index_s = 0; index_s < origSize.u; index_s++) { //horizontal angular
            mFirstPixelPosition = cont * origSize.x*2;
            cont++;

            for (int pixelCount = 0; pixelCount < origSize.x*2; pixelCount++) {

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

                    rgb[icomp][mFirstPixelPosition + pixelCount] = clip(
                            double(pixel[icomp] * clipval), 0.0, (double) clipval);
                }
            }
        }
    }
}


//modificado vertical
//void Prediction::YCbCR2RGBVector(float **yCbCr, const Point4D &origSize, float **rgb, int mPGMScale) {}


//original
void Prediction::write(float **rgb, const Point4D &origSize, int mPGMScale, int start_t, int start_s, const std::string fileName) {
    FILE *mViewFilePointer = fopen(fileName.c_str(), "w");
    if (mViewFilePointer == nullptr) {
        printf("unable to open %s view file for writing\n", fileName.c_str());
        //assert(false);
    }

    int mNumberOfFileBytesPerPixelComponent = (mPGMScale <= 255 ? 1 : 2);

    fprintf(mViewFilePointer, "P6\n%d %d\n%d\n", origSize.x * origSize.u, origSize.y * origSize.v, mPGMScale);

    Point4D it_pos;
    for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {
        for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
            for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {
                for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

                    int pos_out = (it_pos.x) + (it_pos.y * origSize.x) + (it_pos.u * origSize.x * origSize.y)
                                  + (it_pos.v * origSize.x * origSize.y * origSize.u);

                    WritePixelToFile(pos_out, rgb, mPGMScale, mNumberOfFileBytesPerPixelComponent, mViewFilePointer);


                }
            }
        }
    }
}

//write modificado para vertical
void Prediction::writeVector(float **rgb, const Point4D &origSize, int mPGMScale, int start_t, int start_s, const std::string fileName) {
    FILE *mViewFilePointer = fopen(fileName.c_str(), "w");
    if (mViewFilePointer == nullptr) {
        printf("unable to open %s view file for writing\n", fileName.c_str());
        //assert(false);
    }

    int mNumberOfFileBytesPerPixelComponent = (mPGMScale <= 255 ? 1 : 2);

    fprintf(mViewFilePointer,"P6\n%d %d\n%d\n", (origSize.u * origSize.x)*2, origSize.v, mPGMScale); // Horizontal Vector

    Point4D it_pos;

    int pos_out;

    for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
        for (int i = 0; i < 2; ++i) {
            for (it_pos.x = 0; it_pos.x < origSize.x; it_pos.x += 1) {
                for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

                    if (i == 0){
                        pos_out = (it_pos.x) + (it_pos.u * origSize.x)
                                      + (it_pos.v * origSize.x * origSize.u);
                    }
                    else {
                        pos_out = (it_pos.x) + (it_pos.u * origSize.x)
                                  + (it_pos.v * origSize.x * origSize.u) + (origSize.x * origSize.u * origSize.v);
                    }


                    WritePixelToFile(pos_out, rgb, mPGMScale, mNumberOfFileBytesPerPixelComponent, mViewFilePointer);
                }
            }
        }
    }

    /*fprintf(mViewFilePointer,"P6\n%d %d\n%d\n", origSize.u, (origSize.v * origSize.x)*2, mPGMScale); // Vertical vector

    Point4D it_pos;

    int pos_out;

    for (int i = 0; i < 2; ++i) {
        for (it_pos.y = 0; it_pos.y < origSize.y; it_pos.y += 1) {
            for (it_pos.v = 0; it_pos.v < origSize.v; it_pos.v += 1) {
                for (it_pos.u = 0; it_pos.u < origSize.u; it_pos.u += 1) {

                    if (i == 0){
                        pos_out = (it_pos.y) + (it_pos.u * origSize.y)
                                  + (it_pos.v * origSize.y * origSize.u);
                    }
                    else {
                        pos_out = (it_pos.y) + (it_pos.u * origSize.y)
                                  + (it_pos.v * origSize.y * origSize.u) + (origSize.y * origSize.u * origSize.v);
                    }


                    WritePixelToFile(pos_out, rgb, mPGMScale, mNumberOfFileBytesPerPixelComponent, mViewFilePointer);
                }
            }
        }
    }*/
}

void Prediction::WritePixelToFile(int pixelPositionInCache, float **rgb, int mPGMScale, int mNumberOfFileBytesPerPixelComponent, FILE *mViewFilePointer) {

    for (int component_index = 0; component_index < 3; component_index++) {
        int ClippedPixelValue = rgb[component_index][pixelPositionInCache];
        if (ClippedPixelValue > mPGMScale)
            ClippedPixelValue = mPGMScale;
        if (ClippedPixelValue < 0)
            ClippedPixelValue = 0;
        unsigned short bigEndianPixelValue = (mNumberOfFileBytesPerPixelComponent == 2) ? change_endianness_16b(
                ClippedPixelValue) : ClippedPixelValue;


        fwrite(&bigEndianPixelValue, mNumberOfFileBytesPerPixelComponent, 1, mViewFilePointer);
        fflush(mViewFilePointer);
    }

       // std::cout << rgb[2][pixelPositionInCache] << std::endl;

}

unsigned short Prediction::change_endianness_16b(unsigned short val) {
    return (val << 8u) | ((val >> 8u) & 0x00ff);
}

void Prediction::blockGenerator(const Point4D &origSize, int mode, int mPGMScale, int start_t, int start_s, const std::string fileName){
    Point4D it_pos_out;
    float block[origSize.getNSamples()];

    it_pos_out.x = 0;
    it_pos_out.u = 0;
    it_pos_out.v = 0;
    it_pos_out.y = 0;

    float lumaValue = 0;

    if(mode == 1){//Horizontal
        // percorre vetor out na ordem horizontal espacial
        for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {
            for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {

                // percorre vetor out na ordem horizontal angular
                for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {
                    for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {

                        int pos_out =
                                (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                        if((it_pos_out.v == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                           || (it_pos_out.v == 1 && it_pos_out.u <= 1)
                           || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && it_pos_out.u == 0)) { //pixels pretos

                            block[pos_out] = -448;

                        }else{
                            block[pos_out] = lumaValue;
                        }
                    }
                    lumaValue += 20;
                }
                lumaValue = 0;
            }
        }
    }else if(mode == 2){//Vertical
        // percorre vetor out na ordem vertical espacial
        for (it_pos_out.x = 0; it_pos_out.x < origSize.x; it_pos_out.x += 1) {
            for (it_pos_out.y = 0; it_pos_out.y < origSize.y; it_pos_out.y += 1) {

                // percorre vetor out na ordem vertical angular
                for (it_pos_out.u = 0; it_pos_out.u < origSize.u; it_pos_out.u += 1) {
                    for (it_pos_out.v = 0; it_pos_out.v < origSize.v; it_pos_out.v += 1) {

                        int pos_out =
                                (it_pos_out.x) + (it_pos_out.y * origSize.x) + (it_pos_out.u * origSize.x * origSize.y)
                                + (it_pos_out.v * origSize.x * origSize.y * origSize.u);

                        if((it_pos_out.v == 0 && (it_pos_out.u <= 3 || it_pos_out.u > origSize.u - 3))
                           || (it_pos_out.v == 1 && it_pos_out.u <= 1)
                           || ((it_pos_out.v == 2 || it_pos_out.v == 3 || it_pos_out.v > origSize.v - 3) && it_pos_out.u == 0)) { //pixels pretos

                            block[pos_out] = -448;

                        }else{
                            block[pos_out] = lumaValue;
                        }

                    }
                    lumaValue += 20;
                }
                lumaValue = 0;
            }
        }
    }

    float pred[origSize.getNSamples()],
            ref[origSize.getNSamples()];
    float blockVector[(origSize.x * origSize.u * origSize.v)*2];
    //this->angularPrediction(0, 0, block, origSize, pred, 0, ref);
    //this->angularPredictionVector(0, 0, block, origSize, pred, 0, ref, blockVector);

    float *origUp[3], *origRGB[3], *predUp[3], *predRGB[3], *blockVectorUp[3], *blockVectorRGB[3];
    for(int i = 0; i <3; ++i) {
        origUp[i] = new float[origSize.getNSamples()];
        origRGB[i] = new float[origSize.getNSamples()];
        predUp[i] = new float[origSize.getNSamples()];
        predRGB[i] = new float[origSize.getNSamples()];
        blockVectorUp[i] = new float[(origSize.x * origSize.u * origSize.v)*2];
        blockVectorRGB[i] = new float[(origSize.x * origSize.u * origSize.v)*2];
    }
    for(int i = 0; i <3; ++i){
        if(i==0){
            for (int j = 0; j < origSize.getNSamples(); ++j) {
                origUp[i][j] = block[j];
                predUp[i][j] = pred[j];
            }
        }else{
            for (int j = 0; j < origSize.getNSamples(); ++j) {
                origUp[i][j] = 0;
                predUp[i][j] = 0;
            }
        }
    }

    for(int i = 0; i <3; ++i){
        if(i==0){
            for (int j = 0; j < (origSize.x * origSize.u * origSize.v)*2; ++j) {
                blockVectorUp[i][j] = blockVector[j];
            }
        }else{
            for (int j = 0; j < (origSize.x * origSize.u * origSize.v)*2; ++j) {
                blockVectorUp[i][j] = 0;
            }
        }
    }

    this->YCbCR2RGB(origUp, origSize, origRGB, mPGMScale);
    this->YCbCR2RGB(predUp, origSize, predRGB, mPGMScale);

    this->YCbCR2RGBVector(blockVectorUp, origSize, blockVectorRGB, mPGMScale);


    //std::cout << "ind: " << ind << " pos: " << pos << std::endl;

    this->writeVector(blockVectorRGB, origSize, mPGMScale, start_t, start_s, fileName + "_genvector_" + std::to_string(mode) + ".ppm");

    this->write(origRGB, origSize, mPGMScale, start_t, start_s, fileName + "_gen_" + std::to_string(mode) + ".ppm");
    this->write(predRGB, origSize, mPGMScale, start_t, start_s, fileName + "_pred_" + std::to_string(mode) + ".ppm");
}

void Prediction::update(float *curr, bool available, uint blockSize) {
    this->pred_references.erase(this->pred_references.begin());
    this->pred_references.emplace_back(curr, available, blockSize);
}

void Prediction::init_references() {
    for (int i = 0; i < this->resol_x + 1; ++i) this->pred_references.emplace_back();
}

//EDUARDO END