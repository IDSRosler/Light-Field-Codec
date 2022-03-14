#include "SubpartitionModel.h"

SubpartitionModel::SubpartitionModel(const int *bitstream, const Point4D &dim_block, EntropyReport *csv_report) {

    this->block4D = new Block4D(dim_block.x,dim_block.y,dim_block.u,dim_block.v);

    for (int it_y = 0; it_y < this->block4D->dimention.y; ++it_y) {
        for (int it_x = 0; it_x < this->block4D->dimention.x; ++it_x) {
            for (int it_v = 0; it_v < this->block4D->dimention.v; ++it_v) {
                for (int it_u = 0; it_u < this->block4D->dimention.u; ++it_u) {
                    /*if (it_x == 16 || it_y == 16){
                        this->block4D->data[it_x][it_y][it_u][it_v] = 0;
                    }else*/
                        this->block4D->data[it_x][it_y][it_u][it_v] = bitstream[it_x +
                                                                                (it_y * dim_block.x) +
                                                                                (it_u * dim_block.x * dim_block.y) +
                                                                                (it_v * dim_block.x * dim_block.y * dim_block.u)];
                }
            }
        }
    }

    this->root = new NodeBlock({0,0},
                               {(int)this->block4D->dimention.x, (int)this->block4D->dimention.y},
                               this->block4D->dimention, 0);

    this->SetNodeAttributes(this->root, csv_report);

    this->MakeTree(this->root, {0,0}, csv_report);
}

SubpartitionModel::~SubpartitionModel() {
}

void SubpartitionModel::MakeTree(NodeBlock *node, Position middleBefore, EntropyReport *csv_report) {
    if (!node->attributes.has_significant_value || node->dimention.x < 2 || node->dimention.y < 2) {
        return;
    }
    else{
        Position middle = {static_cast<int>(std::ceil(static_cast<float>(node->dimention.x/2.0))),
                           static_cast<int>(std::ceil(static_cast<float>(node->dimention.y/2.0)))};

        for (int i = 0; i < SUBDIVISIONS; ++i) {
            Position pos = this->ComputePositions(i, middleBefore, middle);

            node->child[i] = new NodeBlock(this->startP,
                                           this->endP,
                                           {(uint)(this->endP.x - this->startP.x),
                                            (uint)(this->endP.y - this->startP.y),
                                            node->dimention.u,
                                            node->dimention.v},
                                            node->level + 1);

            this->SetNodeAttributes(node->child[i], csv_report);
            this->MakeTree(node->child[i], pos, csv_report);
        }
    }
}

void SubpartitionModel::SetNodeAttributes(NodeBlock *node, EntropyReport *csv_report) {
    int n_zero = 0,
        n_one = 0,
        n_two = 0,
        n_greater_than_two = 0,
        n_sig_coeff = 0,
        max_value = 0;

    float mean_value = 0.0;

    bool has_significant_value = false;

    for (int it_y = node->start_index.y; it_y < node->end_index.y; ++it_y) {
        for (int it_x = node->start_index.x; it_x < node->end_index.x; ++it_x) {
            for (int it_v = 0; it_v < node->dimention.v; ++it_v) {
                for (int it_u = 0; it_u < node->dimention.u; ++it_u) {
                    if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) > max_value) max_value = abs(this->block4D->data[it_x][it_y][it_u][it_v]);
                    if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) == 0) ++n_zero;
                    else if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) == 1) ++n_one;
                    else if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) == 2) ++n_two;
                    else if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) > 2) ++n_greater_than_two;
                    if (std::abs(this->block4D->data[it_x][it_y][it_u][it_v]) != 0) {++n_sig_coeff; has_significant_value = true;}
                    mean_value += static_cast<float>(abs(this->block4D->data[it_x][it_y][it_u][it_v]));
                }
            }
        }
    }
    mean_value = mean_value / static_cast<float>(node->dimention.getNSamples());
    node->attributes.n_one = n_one;
    node->attributes.n_two = n_two;
    node->attributes.n_zero = n_zero;
    node->attributes.max_value = max_value;
    node->attributes.mean_value = mean_value;
    node->attributes.n_sig_coeff = n_sig_coeff;
    node->attributes.n_greater_than_two = n_greater_than_two;
    node->attributes.has_significant_value = has_significant_value;

    csv_report->writeSubpartitions(
        node->level,
        node->start_index.x,
        node->start_index.y,
        node->end_index.x,
        node->end_index.y,
        node->attributes.n_zero,
        node->attributes.n_one,
        node->attributes.n_two,
        node->attributes.n_greater_than_two,
        node->attributes.max_value,
        node->attributes.mean_value,
        node->attributes.has_significant_value,
        static_cast<float>(node->dimention.getNSamples())
    );
}

Position SubpartitionModel::GetStartPosition(int index, Position middle) {
    if (index == 0) return {0,0};
    if (index == 1) return {middle.x,0};
    if (index == 2) return {0, middle.y};
    if (index == 3) return {middle.x,middle.y};
}

Position SubpartitionModel::ComputePositions(int index, Position middleBefore, Position middle) {
    Position pos = this->GetStartPosition(index, middle);

    this->startP.x = pos.x + middleBefore.x;
    this->startP.y = pos.y + middleBefore.y;

    this->endP.x = this->startP.x + middle.x;
    this->endP.y = this->startP.y + middle.y;

    if(this->endP.x > this->root->dimention.x) this->endP.x = this->root->dimention.x;
    if(this->endP.y > this->root->dimention.y) this->endP.y = this->root->dimention.y;

    return pos;
}

void SubpartitionModel::_deleteTree(NodeBlock *node) {
    if (node == nullptr) return;

    for (int i = 0; i < node->child.size(); ++i) {
        this->_deleteTree(node->child[i]);
        node->child.clear();
    }
    delete node;
}

void SubpartitionModel::DeleteTree() {
    delete [] this->block4D->data;

    this->_deleteTree(this->root);
    this->root = nullptr;
}