#ifndef LF_CODEC_TREE_H
#define LF_CODEC_TREE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <queue>
#include <fstream>
#include <cassert>

#include "EncoderParameters.h"
#include "Point4D.h"
#include "Typedef.h"

#define HEXADECA 16

using namespace std;

struct SyntacticElements {
    int last;                   // last significant element
    vector<short int> sig;      // significant elements after the last
    vector<short int> gr_one;   // indicates whether significant elements have an absolute value greater than 1
    vector<short int> gr_two;   // indicates whether significant elements have an absolute value greater than 2
    vector<short int> sign;     // sign of  significant elements (1 - negative | 0 - positive)
    vector<int> rem;            // rest of the significant elements greater than one (and two)

    void reset() {
        this->last = -1;
        this->sig.clear();
        this->gr_one.clear();
        this->gr_two.clear();
        this->sign.clear();
        this->rem.clear();
    }
};

struct Point_4D{
    int x = 0, y = 0, u = 0, v = 0;
};

struct Hypercube{
    int ****data;
    Point_4D dim = {0,0,0,0};

    Hypercube(int x, int y, int u, int v){
        this->dim.x = x;
        this->dim.y = y;
        this->dim.u = u;
        this->dim.v = v;

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

struct Attributes{
    uint n_zero = 0,
         n_one = 0,
         n_two = 0,
         n_greater_than_two = 0,
         hypercubo_size = 0,
         sig_coeff = 0;


    int max_value = 0;

    float mean_value = 0.0;

    bool significant_value = false;
};

struct Node{
    Point_4D start{0,0,0,0};
    Point_4D end{0,0,0,0};
    Point_4D hypercube_dim{0,0,0,0};
    Point_4D node_pos{0,0,0,0};

    Attributes *att;

    vector<Node *> child;

    Node(Point_4D start, Point_4D end, Point_4D hypercube_dim){
        this->start = start;
        this->end = end;
        this->hypercube_dim = hypercube_dim;

        this->att = nullptr;

        for (int i = 0; i < HEXADECA; ++i) {
            this->child.push_back(nullptr);
        }
    }

    void SetNodePosition(Point_4D pos){
        this->node_pos = pos;
    }

    void SetAttributes(Attributes *attr){
        this->att = attr;
    }
};

class Tree {
public:
    Tree();
    Node* CreateRoot(const int *bitstream, const Point4D &dim_block);

    void ComputeLast(int &last);
    void ComputeRun(vector<int> &v_run, int last);
    void ComputeSyntacticElements(vector<SyntacticElements> &lfbpu_elements, int last);

    void CreateTree(Node * root, const Point4D &pos, Point_4D middle_before);
    void DeleteTree(Node** node_ref);

    ~Tree();

private:
    void ComputeAttributes(Node *node, int start_x, int end_x, int start_y, int end_y, int start_u, int end_u, int start_v, int end_v);
    void ComputePositions(Point_4D start, Point_4D middle_before, Point_4D middle);
    Point_4D ComputeStart(int index, Point_4D middle);

    void HypercubePosition(Point_4D *middle);
    void _deleteTree(Node* node);

    void SortBufferPositions();
    void LFBPUToVector(vector<int> &v_coefficients, int index);

    int index_sorted[256] = {
            0,1,4,5,16,17,20,21,64,65,68,69,80,81,84,85,2,3,6,7,18,19,22,23,66,67,70,71,82,83,86,87,8,9,12,13,24,25,28,29,72,73,76,77,88,89,92,93,10,
            11,14,15,26,27,30,31,74,75,78,79,90,91,94,95,32,33,36,37,48,49,52,53,96,97,100,101,112,113,116,117,34,35,38,39,50,51,54,55,98,99,102,103,
            114,115,118,119,40,41,44,45,56,57,60,61,104,105,108,109,120,121,124,125,42,43,46,47,58,59,62,63,106,107,110,111,122,123,126,127,128,129,
            132,133,144,145,148,149,192,193,196,197,208,209,212,213,130,131,134,135,146,147,150,151,194,195,198,199,210,211,214,215,136,137,140,141,
            152,153,156,157,200,201,204,205,216,217,220,221,138,139,142,143,154,155,158,159,202,203,206,207,218,219,222,223,160,161,164,165,176,177,
            180,181,224,225,228,229,240,241,244,245,162,163,166,167,178,179,182,183,226,227,230,231,242,243,246,247,168,169,172,173,184,185,188,189,
            232,233,236,237,248,249,252,253,170,171,174,175,186,187,190,191,234,235,238,239,250,251,254,255
    };

    vector<Node *> order4_SubPartitionsBuffer;

    Hypercube *hypercube = nullptr;

    Point_4D next_start_position = {0,0,0,0};
    Point_4D next_end_position = {0,0,0,0};
    Point_4D hy_pos = {0,0,0,0};

    uint size = 0;
};

#endif //LF_CODEC_TREE_H
