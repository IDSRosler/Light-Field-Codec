//
// Created by igor on 09/04/2022.
//

#ifndef LF_CODEC_ENTROPYSTRUCTSMODEL_H
#define LF_CODEC_ENTROPYSTRUCTSMODEL_H

#define SUBDIVISIONS 4

#include <queue>

#include "Point4D.h"

struct Position{
    int x, y;
};

struct Block4D{
    Point4D dimention = {0,0,0,0};
    int ****data;

    Block4D(int x, int y, int u, int v){
        this->dimention.x = x;
        this->dimention.y = y;
        this->dimention.u = u;
        this->dimention.v = v;

        this->data = (int ****) calloc(x, sizeof(int ***));
        for (int i = 0; i < x; ++i) {
            this->data[i] = (int ***) calloc(y, sizeof(int **));
            for (int j = 0; j < y; ++j) {
                this->data[i][j] = (int **) calloc(u, sizeof(int *));
                for (int k = 0; k < u; ++k) {
                    this->data[i][j][k] = (int *) calloc(v, sizeof(int));
                }
            }
        }
    }
};

struct NodeAttributes{
    int n_zero = 0,
            n_one = 0,
            n_two = 0,
            n_greater_than_two = 0,
            n_sig_coeff = 0,
            max_value = 0;

    float mean_value = 0.0;

    bool has_significant_value = false;
};

struct NodeBlock{
    Position start_index{0,0};
    Position end_index{0,0};
    Point4D dimention{0,0,0,0};
    int level;

    std::vector<NodeBlock *> child;

    NodeAttributes attributes;

    NodeBlock(Position start, Position end, Point4D dim, int level){
        this->start_index = start;
        this->end_index = end;
        this->dimention = dim;
        this->level = level;

        for (int i = 0; i < SUBDIVISIONS; ++i) {
            this->child.push_back(nullptr);
        }
    }
};

struct LRE_Struct{
    int level;
    int run;
};

struct Syntactic_Elements{
    int last_sig_coeff_u=-1;
    int last_sig_coeff_v=-1;
    std::queue<bool> sig_coeff_flag;
    std::queue<bool> coeff_abs_level_greater1_flag;
    std::queue<bool> coeff_abs_level_greater2_flag;
    std::queue<uint> coeff_abs_level_remaining;
    std::queue<bool> coeff_sign_flag;
};

#endif // LF_CODEC_ENTROPYSTRUCTSMODEL_H