#ifndef LF_CODEC_SUBPARTITIONMODEL_H
#define LF_CODEC_SUBPARTITIONMODEL_H

#include <vector>
#include <cmath>

#include "EntropyReport.h"
#include "ArithmeticStructures.h"

class SubpartitionModel {
  public:
    SubpartitionModel(const int *bitstream, const Point4D &dim_block, std::vector<bool> &treeFlags, EntropyReport *csv_report, int hy);
    ~SubpartitionModel();
    void DeleteTree();

private:
    void _deleteTree(NodeBlock *node);
    void SetNodeAttributes(NodeBlock *node, EntropyReport *csv_report);
    void MakeTree(NodeBlock *node, std::vector<bool> &treeFlags, EntropyReport *csv_report, int hy);
    Position GetStartPosition(int index, NodeBlock *node, Position middle);
    Position GetEndPosition(int index, NodeBlock *node, Position middle);

    Block4D *block4D = nullptr;
    NodeBlock *root = nullptr;
    Position startP = {0,0};
    Position endP = {0,0};
};

#endif // LF_CODEC_SUBPARTITIONMODEL_H
