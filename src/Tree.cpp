#include "Tree.h"

Tree::Tree() = default;

Tree::~Tree() = default;

Node* Tree::CreateRoot(const int *bitstream, const Point4D &dim_block) {

    //this->hypercube = new Hypercube(dim_block.x, dim_block.y, dim_block.u, dim_block.v); //Not zero padding

    this->hypercube = new Hypercube(16,16,16,16);

    this->order4_SubPartitionsBuffer.clear();

    int i=0;
    for (int it_v = 0; it_v < this->hypercube->dim.v; ++it_v) {
        for (int it_u = 0; it_u < this->hypercube->dim.u; ++it_u) {
            for (int it_y = 0; it_y < this->hypercube->dim.y; ++it_y) {
                for (int it_x = 0; it_x < this->hypercube->dim.x; ++it_x) {
                    if (i < static_cast<int>(dim_block.x * dim_block.y * dim_block.u * dim_block.v)) {
                        this->hypercube->data[it_x][it_y][it_u][it_v] = bitstream[i];
                        i++;
                    }
                    else this->hypercube->data[it_x][it_y][it_u][it_v] = 0;
                }
            }
        }
    }

    Node *root = new Node({0,0,0,0}, this->hypercube->dim, this->hypercube->dim); //Zero padding

    return root;
}

void Tree::CreateTree(Node * root, const Point4D &pos, Point_4D middle_before) {
    if (HEXADECA_TREE_PARTITION == 0)
        this->size = 16;
    else if (HEXADECA_TREE_PARTITION == 1)
        this->size = 8;
    else if (HEXADECA_TREE_PARTITION == 2)
        this->size = 4;
    else return;

    if (root->hypercube_dim.x <= this->size || root->hypercube_dim.y <= this->size || root->hypercube_dim.u <= this->size || root->hypercube_dim.v <= this->size) {

        this->ComputeAttributes(root, root->start.x, root->end.x, root->start.y, root->end.y, root->start.u, root->end.u, root->start.v, root->end.v);

        root->SetNodePosition(this->hy_pos);

        this->order4_SubPartitionsBuffer.push_back(root);

        this->hy_pos = {0,0,0,0};

        return;
    }
    else{
        this->ComputeAttributes(root, root->start.x, root->end.x, root->start.y, root->end.y, root->start.u, root->end.u, root->start.v, root->end.v);

        Point_4D middle = {(int)ceil((double)root->hypercube_dim.x/2),
                           (int)ceil((double)root->hypercube_dim.y/2),
                           (int)ceil((double)root->hypercube_dim.u/2),
                           (int)ceil((double)root->hypercube_dim.v/2)};

        Point_4D start;

        for (int i = 0; i < HEXADECA; ++i) {
            start = this->ComputeStart(i, middle);

            this->ComputePositions(start, middle_before, middle);
            this->HypercubePosition(&middle);

            root->child[i] = new Node(this->next_start_position, this->next_end_position, middle);
            this->CreateTree(root->child[i], pos, start);
        }
    }
}

Point_4D Tree::ComputeStart(int index, Point_4D middle) {
    if (index == 0) return {0,0,0,0};
    if (index == 1) return {middle.x,0,0,0};
    if (index == 2) return {0,middle.y,0,0};
    if (index == 3) return {middle.x,middle.y,0,0};
    if (index == 4) return {0,0,middle.u,0};
    if (index == 5) return {middle.x,0,middle.u,0};
    if (index == 6) return {0,middle.y,middle.u,0};
    if (index == 7) return {middle.v,middle.y,middle.u,0};
    if (index == 8) return {0,0,0,middle.v};
    if (index == 9) return {middle.x,0,0,middle.v};
    if (index == 10) return {0,middle.y,0,middle.v};
    if (index == 11) return {middle.x,middle.y,0,middle.v};
    if (index == 12) return {0,0,middle.u,middle.v};
    if (index == 13) return {middle.x,0,middle.u,middle.v};
    if (index == 14) return {0,middle.y,middle.u,middle.v};
    if (index == 15) return {middle.x,middle.y,middle.u,middle.v};

    return {0,0,0,0,};
}

void Tree::ComputePositions(Point_4D start, Point_4D middle_before, Point_4D middle) {
    this->next_start_position.x = start.x + middle_before.x;
    this->next_start_position.y = start.y + middle_before.y;
    this->next_start_position.u = start.u + middle_before.u;
    this->next_start_position.v = start.v + middle_before.v;

    this->next_end_position.x = this->next_start_position.x + middle.x;
    this->next_end_position.y = this->next_start_position.y + middle.y;
    this->next_end_position.u = this->next_start_position.u + middle.u;
    this->next_end_position.v = this->next_start_position.v + middle.v;
}

void Tree::ComputeAttributes(Node* node, int start_x, int end_x, int start_y, int end_y, int start_u, int end_u, int start_v, int end_v) {
    float acc = 0.0;
    auto *att = new Attributes();
    for (int it_v = start_v; it_v < end_v; ++it_v) {
        for (int it_u = start_u; it_u < end_u; ++it_u) {
            for (int it_y = start_y; it_y < end_y; ++it_y) {
                for (int it_x = start_x; it_x < end_x; ++it_x) {
                    if (abs(att->max_value) < abs(this->hypercube->data[it_x][it_y][it_u][it_v])) { att->max_value = abs(this->hypercube->data[it_x][it_y][it_u][it_v]);}
                    if (this->hypercube->data[it_x][it_y][it_u][it_v] != 0) {att->significant_value = true; ++att->sig_coeff;}
                    if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 0) {++att->n_zero;}
                    else if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 1) {++att->n_one;}
                    else if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 2) {++att->n_two;}
                    else {++att->n_greater_than_two;}
                    acc += static_cast<float>(abs(this->hypercube->data[it_x][it_y][it_u][it_v]));
                }
            }
        }
    }
    att->hypercubo_size = node->hypercube_dim.x * node->hypercube_dim.y * node->hypercube_dim.u * node->hypercube_dim.v;
    att->mean_value = acc / static_cast<float>(att->hypercubo_size);
    node->SetAttributes(att);
}

void Tree::HypercubePosition(Point_4D *middle) {
    this->hy_pos.x = (this->next_start_position.x != 0) ? this->next_start_position.x / middle->x : this->next_start_position.x,
    this->hy_pos.y = (this->next_start_position.y != 0) ? this->next_start_position.y / middle->y : this->next_start_position.y,
    this->hy_pos.u = (this->next_start_position.u != 0) ? this->next_start_position.u / middle->u : this->next_start_position.u,
    this->hy_pos.v = (this->next_start_position.v != 0) ? this->next_start_position.v / middle->v : this->next_start_position.v;
}

void Tree::ComputeLast(int &last) {
    int index;

    this->SortBufferPositions();

    for (index = static_cast<int>(this->order4_SubPartitionsBuffer.size() - 1); index >= 0; --index) {
        if (this->order4_SubPartitionsBuffer[index]->att->significant_value){
            break;
        }
    }
    if (index < 0){
        index = 0;
    }
    last = index;
}

void Tree::subpartitionReport(EntropyReport &report) {
    for (int index = static_cast<int>(this->order4_SubPartitionsBuffer.size() - 1); index >= 0; --index) {
        report.writeSubpartitions(
          index,
          this->order4_SubPartitionsBuffer[index]->node_pos.x,
          this->order4_SubPartitionsBuffer[index]->node_pos.y,
          this->order4_SubPartitionsBuffer[index]->node_pos.u,
          this->order4_SubPartitionsBuffer[index]->node_pos.v,
          this->order4_SubPartitionsBuffer[index]->att->hypercubo_size,
          this->order4_SubPartitionsBuffer[index]->att->n_zero,
          this->order4_SubPartitionsBuffer[index]->att->n_one,
          this->order4_SubPartitionsBuffer[index]->att->n_two,
          this->order4_SubPartitionsBuffer[index]->att->n_greater_than_two,
          this->order4_SubPartitionsBuffer[index]->att->max_value,
          this->order4_SubPartitionsBuffer[index]->att->mean_value,
          this->order4_SubPartitionsBuffer[index]->att->significant_value
        );
    }
}

void Tree::ComputeRun(vector<int> &v_run, int last) {
    int run = 0;
    int index = last-1;

    while (index >= 0){ // Calcula a corrida de blocos zerados
        while (index >= 0 && !this->order4_SubPartitionsBuffer[index]->att->significant_value){
            ++run;
            --index;
        }
        if (index >= 0){
            --index;
        }
        v_run.push_back(run);
        run = 0;
    }
}

void Tree::ComputeSyntacticElements(vector<SyntacticElements> &lfbpu_elements, int last) {
    int j;
    vector<int> v_coefficients;
    SyntacticElements elem;
    elem.reset();
    v_coefficients.clear();
    for (int i = last; i >= 0; --i) {
        this->LFBPUToVector(v_coefficients, i);
        if (!v_coefficients.empty()){
            for (j = static_cast<int>(v_coefficients.size() - 1); j >= 0; --j) { // compute last (coefficient level)
                if (v_coefficients[j] != 0) {
                    break;
                }
            }
            elem.last = j;
            for (int k = elem.last; k >= 0; --k) { // compute sig (coefficient level)
                if (v_coefficients[k] == 0){
                    elem.sig.push_back(0);
                }else {
                    elem.sig.push_back(1);
                    if (abs(v_coefficients[k]) > 1){ // > 1
                        elem.gr_one.push_back(1);
                        if (abs(v_coefficients[k]) > 2){ // > 2
                            elem.gr_two.push_back(1);
                            elem.rem.push_back(abs(v_coefficients[k]) - 3); //rem
                        } else { // = 2
                            elem.gr_two.push_back(0);
                        }
                    }else { // = 1
                        elem.gr_one.push_back(0);
                    }
                    if (v_coefficients[k] < 0) { // signal
                        elem.sign.push_back(1);
                    } else {
                        elem.sign.push_back(0);
                    }
                }
            }
            v_coefficients.clear();
            lfbpu_elements.push_back(elem);
            elem.reset();
        }
    }
}

void Tree::LFBPUToVector(vector<int> &v_coefficients, int index) {
    if(order4_SubPartitionsBuffer[index]->att->significant_value) {
        for (int v = order4_SubPartitionsBuffer[index]->start.v ; v < order4_SubPartitionsBuffer[index]->end.v ; ++v) {
            for (int u = order4_SubPartitionsBuffer[index]->start.u ; u < order4_SubPartitionsBuffer[index]->end.u ; ++u) {
                for (int y = order4_SubPartitionsBuffer[index]->start.y ; y < order4_SubPartitionsBuffer[index]->end.y ; ++y) {
                    for (int x = order4_SubPartitionsBuffer[index]->start.x ; x < order4_SubPartitionsBuffer[index]->end.x ; ++x) {
                        v_coefficients.push_back(this->hypercube->data[x][y][u][v]);
                    }
                }
            }
        }
    }
}

void Tree::SortBufferPositions() {
    vector<Node *> temp;
    for (int i = 0; i < this->order4_SubPartitionsBuffer.size(); ++i) {
        temp.push_back(nullptr);
    }
    for (int i = 0; i < this->order4_SubPartitionsBuffer.size(); ++i) {
        temp[this->index_sorted[i]] = this->order4_SubPartitionsBuffer[i];
    }
    this->order4_SubPartitionsBuffer = temp;
    temp.clear();
}

void Tree::_deleteTree(Node* node)
{
    if (node == nullptr) return;

    for (int i = 0; i < node->child.size(); ++i) {
        _deleteTree(node->child[i]);
        node->child.clear();
    }
    delete node;
}

void Tree::DeleteTree(Node** node_ref)
{
    delete [] this->hypercube->data;

    _deleteTree(*node_ref);
    *node_ref = nullptr;

    this->order4_SubPartitionsBuffer.clear();
}