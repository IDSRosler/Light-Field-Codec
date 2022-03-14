#ifndef LF_CODEC_SUBPARTITIONMODEL_H
#define LF_CODEC_SUBPARTITIONMODEL_H

#include <vector>
#include <cmath>

#include "Point4D.h"
#include "EntropyReport.h"

#define SUBDIVISIONS 4

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

class SubpartitionModel {
  public:
    SubpartitionModel(const int *bitstream, const Point4D &dim_block, EntropyReport *csv_report);
    ~SubpartitionModel();
    void DeleteTree();

  private:
    void MakeTree(NodeBlock *node, Position middleBefore, EntropyReport *csv_report);
    void _deleteTree(NodeBlock *node);
    void SetNodeAttributes(NodeBlock *node, EntropyReport *csv_report);

    Position ComputePositions(int index, Position middleBefore, Position middle);
    Position GetStartPosition(int index, Position middle);

    Block4D *block4D = nullptr;
    NodeBlock *root = nullptr;
    Position startP = {0,0};
    Position endP = {0,0};
};

#endif // LF_CODEC_SUBPARTITIONMODEL_H
