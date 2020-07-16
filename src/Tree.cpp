#include "Tree.h"

Tree::Tree() {
}

Tree::~Tree() {
}

Node* Tree::CreateRoot(string light_field, uint hypercubo, uint channel, int *bitstream, const Point4D &dim_block) {
    this->props = new LF_Props(light_field, hypercubo, channel);

    //this->hypercube = new Hypercube(dim_block.x, dim_block.y, dim_block.u, dim_block.v); //Not zero padding

    this->hypercube = new Hypercube(16,16,16,16);

    this->order4_SubPartitionsBuffer.clear();
    this->order8_SubPartitionsBuffer.clear();

    this->id = 0;
    this->CBF_bits_per_hypercube = 0;

    int i=0;
    for (int it_v = 0; it_v < this->hypercube->dim.v; ++it_v) {
        for (int it_u = 0; it_u < this->hypercube->dim.u; ++it_u) {
            for (int it_y = 0; it_y < this->hypercube->dim.y; ++it_y) {
                for (int it_x = 0; it_x < this->hypercube->dim.x; ++it_x) {
                    if (i < (dim_block.x * dim_block.y * dim_block.u * dim_block.v)) {
                        this->hypercube->data[it_x][it_y][it_u][it_v] = bitstream[i];
                        i++;
                    }
                    else this->hypercube->data[it_x][it_y][it_u][it_v] = 0;
                }
            }
        }
    }

    Node *root = new Node({0,0,0,0}, this->hypercube->dim, this->hypercube->dim, 0, -1); //Zero padding

    return root;
}

void Tree::CreateTree(Node * root, uint level, const Point4D &pos, const Point_4D &hypercubo_pos, Point_4D middle_before) {
    if (HEXADECA_TREE_PARTITION == 0)
        this->size = 16;
    else if (HEXADECA_TREE_PARTITION == 1)
        this->size = 8;
    else if (HEXADECA_TREE_PARTITION == 2)
        this->size = 4;
    else return;

    if (root->hypercube_dim.x <= this->size || root->hypercube_dim.y <= this->size || root->hypercube_dim.u <= this->size || root->hypercube_dim.v <= this->size) { //order 4 se HEXADECA_TREE_PARTITION == 2
        if (HEXADECA_TREE_PARTITION == 0){
            this->next_start_position.x = pos.x;
            this->next_start_position.y = pos.y;
            this->next_start_position.u = pos.u;
            this->next_start_position.v = pos.v;
        }
        else if (HEXADECA_TREE_PARTITION == 1 || HEXADECA_TREE_PARTITION == 2){
            this->next_start_position.x += pos.x;
            this->next_start_position.y += pos.y;
            this->next_start_position.u += pos.u;
            this->next_start_position.v += pos.v;
        }
#if HEXADECA_TREE_TYPE == 0
        this->ComputeAttributes(root, root->start.x, root->end.x, root->start.y, root->end.y, root->start.u, root->end.u, root->start.v, root->end.v);

#if HEXADECA_TREE_CODEC_MODE == 0
        if (root->att->significant_value == 1 || root->att->significant_value == 0 ){
            ++this->CBF_bits_per_hypercube;
        }
#endif

        root->SetNodePosition(this->hy_pos);

#if HEXADECA_TREE_CODEC_MODE == 1 || HEXADECA_TREE_CODEC_MODE == 2 || HEXADECA_TREE_CODEC_MODE == 4
        this->order4_SubPartitionsBuffer.push_back(root);
#endif

#if HEXADECA_TREE_CODEC_MODE == 5
        this->WriteAttributesInFile(level, this->hy_pos, root);
#endif

        this->hy_pos = {0,0,0,0};
#endif
#if HEXADECA_TREE_TYPE == 1
        if(this->props->hypercubo == HYPERCUBE_POSITION && this->props->channel == HYPERCUBE_CHANNEL)
            this->ComputeValues(root, root->start.x, root->end.x, root->start.y, root->end.y, root->start.u, root->end.u, root->start.v, root->end.v, level, this->hy_pos);
        this->hy_pos = {0,0,0,0};
#endif
        return;
    }
    else{
#if HEXADECA_TREE_TYPE == 0
        this->ComputeAttributes(root, root->start.x, root->end.x, root->start.y, root->end.y, root->start.u, root->end.u, root->start.v, root->end.v);

#if HEXADECA_TREE_CODEC_MODE == 3
        //this->WriteAttributesInFile(level, this->hy_pos, root); //prox_pos
#endif

#endif
        Point_4D middle = {(int)ceil((double)root->hypercube_dim.x/2),
                           (int)ceil((double)root->hypercube_dim.y/2),
                           (int)ceil((double)root->hypercube_dim.u/2),
                           (int)ceil((double)root->hypercube_dim.v/2)};

        Point_4D start = {0,0,0,0};

        uint next_level = level + 1;

#if HEXADECA_TREE_CODEC_MODE == 3
        if(level == 1) {
            this->order8_SubPartitionsBuffer.push_back(root);
        }
#endif

#if HEXADECA_TREE_CODEC_MODE == 0
        if (root->att->significant_value == 1){ // original ou ordem 8
            ++this->CBF_bits_per_hypercube;
        }
        else if(root->att->significant_value == 0){
            ++this->CBF_bits_per_hypercube;
            return;
        }
#endif

        for (int i = 0; i < HEXADECA; ++i) {
            this->id++;
            start = this->ComputeStart(i, middle);

            this->ComputePositions(start, middle_before, middle);
            this->HypercubePosition(&middle);

            root->child[i] = new Node(this->next_start_position, this->next_end_position, middle, this->id, root->id);
            this->CreateTree(root->child[i], next_level, pos, hypercubo_pos, start);
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
    int acc = 0;
    Attributes *att = new Attributes();
    for (int it_v = start_v; it_v < end_v; ++it_v) {
        for (int it_u = start_u; it_u < end_u; ++it_u) {
            for (int it_y = start_y; it_y < end_y; ++it_y) {
                for (int it_x = start_x; it_x < end_x; ++it_x) {
                    acc += abs(this->hypercube->data[it_x][it_y][it_u][it_v]);
                    if (abs(att->max_value) < abs(this->hypercube->data[it_x][it_y][it_u][it_v])) { att->max_value = abs(this->hypercube->data[it_x][it_y][it_u][it_v]);}
                    if (this->hypercube->data[it_x][it_y][it_u][it_v] != 0) {att->significant_value = true;}
                    if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 0) {++att->n_zero;}
                    else if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 1) {++att->n_one;}
                    else if (abs(this->hypercube->data[it_x][it_y][it_u][it_v]) == 2) {++att->n_two;}
                    else {++att->n_greater_than_two;}
                }
            }
        }
    }
    att->hypercubo_size = node->hypercube_dim.x * node->hypercube_dim.y * node->hypercube_dim.u * node->hypercube_dim.v;
    att->mean_value = (float)acc / att->hypercubo_size;
    node->SetAttributes(att);
/*
    if(att->significant_value){
        ++this->CBF_bits_per_hypercube;
    }*/
}

void Tree::ComputeValues(Node* node, int start_x, int end_x, int start_y, int end_y, int start_u, int end_u, int start_v, int end_v, uint level, Point_4D &pos) {
    Point_4D position = {0,0,0,0};
    for (int it_v = start_v; it_v < end_v; ++it_v) {
        for (int it_u = start_u; it_u < end_u; ++it_u) {
            for (int it_y = start_y; it_y < end_y; ++it_y) {
                for (int it_x = start_x; it_x < end_x; ++it_x) {
                    position = {it_x, it_y, it_u, it_v};
                    this->WriteValuesInFile(level, position, position, abs(this->hypercube->data[it_x][it_y][it_u][it_v]));
                }
            }
        }
    }
}

void Tree::HypercubePosition(Point_4D *middle) {
    this->hy_pos.x = (this->next_start_position.x != 0) ? this->next_start_position.x / middle->x : this->next_start_position.x,
    this->hy_pos.y = (this->next_start_position.y != 0) ? this->next_start_position.y / middle->y : this->next_start_position.y,
    this->hy_pos.u = (this->next_start_position.u != 0) ? this->next_start_position.u / middle->u : this->next_start_position.u,
    this->hy_pos.v = (this->next_start_position.v != 0) ? this->next_start_position.v / middle->v : this->next_start_position.v;
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

/* Deletes a tree and sets the root as NULL */
void Tree::DeleteTree(Node** node_ref)
{
    delete [] this->hypercube->data;

    _deleteTree(*node_ref);
    *node_ref = nullptr;
}

void Tree::OpenFile(string path) {
    this->file.open(path + "HexadecaTree.csv");

    assert(this->file.is_open());

    this->SetFileAttributs();
}

void Tree::SetFileAttributs(){
#if HEXADECA_TREE_CODEC_MODE == 0 //CBF hierarquico
        this->file <<
           "Light_Field" << SEP <<
           "Hypercubo" << SEP <<
           "Channel" << SEP <<
           "CBF_Bits_Per_Hypercube" << SEP << endl;
#elif HEXADECA_TREE_CODEC_MODE == 1 //LAST + CBF
        this->file <<
           "Light_Field" << SEP <<
           "Hypercubo" << SEP <<
           "Channel" << SEP <<
           "Last_CBF" << SEP << endl;
#elif HEXADECA_TREE_CODEC_MODE == 2 || HEXADECA_TREE_CODEC_MODE == 4  //LAST + RUN
        this->file <<
           "Light_Field" << SEP <<
           "Hypercubo" << SEP <<
           "Channel" << SEP <<
           "Last_Run" << SEP << endl;
#elif HEXADECA_TREE_CODEC_MODE == 3  //Multiple Levels of LAST + RUN
    this->file <<
           "Light_Field" << SEP <<
           "Hypercubo" << SEP <<
           "Channel" << SEP <<
           "Multiple_Levels_Of_Last_Run" << SEP << endl;
#else //base
        this->file <<
           "Light_Field" << SEP <<
#if HEXADECA_TREE_TYPE == 0
           "Parent" << SEP <<
           "Id" << SEP <<
#endif
           "Partition" << SEP <<
           "Hypercubo" << SEP <<
           "Channel" << SEP <<
           "Pos_x" << SEP <<
           "Pos_y" << SEP <<
           "Pos_u" << SEP <<
           "Pos_v" << SEP <<
#if HEXADECA_TREE_TYPE == 0
           "Hypercubo_Size" << SEP <<
           "N_Zero" << SEP <<
           "N_One" << SEP <<
           "N_Two" << SEP <<
           "N_Greater_Than_Two" << SEP <<
           "Abs_Max_Value" << SEP <<
           "Abs_Mean_value" << SEP <<
           "Significant_Value" << SEP << endl;
#elif HEXADECA_TREE_TYPE == 1
           "X" << SEP <<
           "Y" << SEP <<
           "U" << SEP <<
           "V" << SEP <<
           "Valor" << SEP << endl;
#endif
#endif
}

void Tree::CloseFile() {
    this->file.close();
}

void Tree::WriteAttributesInFile(uint level, Point_4D &pos, Node* node){
    string order = "";
    order = (level == 0) ? "Original" : (level == 1) ? "Order_8" : "Order_4";
    this->file <<
         this->props->light_field_name << SEP <<
         node->parent << SEP <<
         node->id << SEP <<
         order << SEP <<
         this->props->hypercubo << SEP <<
         this->props->channel << SEP <<
         pos.x << SEP <<
         pos.y << SEP <<
         pos.u << SEP <<
         pos.v << SEP <<
         node->att->hypercubo_size << SEP <<
         node->att->n_zero << SEP <<
         node->att->n_one << SEP <<
         node->att->n_two << SEP <<
         node->att->n_greater_than_two <<  SEP <<
         node->att->max_value << SEP <<
         node->att->mean_value << SEP <<
         node->att->significant_value << SEP << endl;
}

void Tree::WriteValuesInFile(uint level, Point_4D &pos, Point_4D &position, int value){
    string order = "";
    order = (level == 0) ? "Original" : (level == 1) ? "Order_8" : "Order_4";
    this->file <<
         this->props->light_field_name << SEP <<
         order << SEP <<
         this->props->hypercubo << SEP <<
         this->props->channel << SEP <<
         pos.x << SEP <<
         pos.y << SEP <<
         pos.u << SEP <<
         pos.v << SEP <<
         position.x << SEP <<
         position.y << SEP <<
         position.u << SEP <<
         position.v << SEP <<
         value << SEP << endl;
}


void Tree::ComputeHierarchicalCBF(){
    this->file <<
        this->props->light_field_name << SEP <<
        this->props->hypercubo << SEP <<
        this->props->channel << SEP <<
        this->CBF_bits_per_hypercube << SEP << endl;
}

void Tree::ComputeLastCBF() {
    this->SortBufferPositions();
    while (this->order4_SubPartitionsBuffer.size() > 0){
        if (this->order4_SubPartitionsBuffer.back()->att->significant_value != 0){
            /*cout << "Last sub-hypercube with significant value: (x: " << this->order4_SubPartitionsBuffer.back()->node_pos.x <<
                    ",y: " << this->order4_SubPartitionsBuffer.back()->node_pos.y <<
                    ",u: " << this->order4_SubPartitionsBuffer.back()->node_pos.u <<
                    ",v: " << this->order4_SubPartitionsBuffer.back()->node_pos.v << ")" << endl;*/
            break;
        }else{
            this->order4_SubPartitionsBuffer.pop_back();
        }
    }
    this->file <<
               this->props->light_field_name << SEP <<
               this->props->hypercubo << SEP <<
               this->props->channel << SEP <<
               (8 + this->order4_SubPartitionsBuffer.size()) << SEP << endl; // 8 bits do Last + 1 flag para cada subpartição restante
}

void Tree::ComputeLastRun() {
    int run = 0;
    int bits_run = 0;
    this->SortBufferPositions();
    while (this->order4_SubPartitionsBuffer.size() > 0){ // Encontra o Last no buffer de ordem 4
        if (this->order4_SubPartitionsBuffer.back()->att->significant_value != 0){
            break;
        }else{
            this->order4_SubPartitionsBuffer.pop_back();
        }
    }
    while (this->order4_SubPartitionsBuffer.size() > 0){ // Calcula a corrida de blocos zerados
        while(this->order4_SubPartitionsBuffer.size() > 0 && this->order4_SubPartitionsBuffer.back()->att->significant_value == 0){
            run++;
            this->order4_SubPartitionsBuffer.pop_back();
        }
        if (this->order4_SubPartitionsBuffer.size() > 0) {
            this->order4_SubPartitionsBuffer.pop_back();
        }
        bits_run += this->BitsExpGolomb(run);
        run = 0;
    }
    this->file <<
               this->props->light_field_name << SEP <<
               this->props->hypercubo << SEP <<
               this->props->channel << SEP <<
               (8 + bits_run) << SEP << endl; // 8 bits do Last + bits referente a corrida

}

void Tree::ComputeLastRun_0_1() {
    int run_0 = 0;
    int run_1 = 0;
    int bits_run = 0;
    this->SortBufferPositions();
    while (this->order4_SubPartitionsBuffer.size() > 0){ // Encontra o Last no buffer de ordem 4
        if (this->order4_SubPartitionsBuffer.back()->att->significant_value != 0){
            break;
        }else{
            this->order4_SubPartitionsBuffer.pop_back();
        }
    }
    while (this->order4_SubPartitionsBuffer.size() > 0){ // Calcula a corrida de blocos zerados
        while(this->order4_SubPartitionsBuffer.size() > 0 && this->order4_SubPartitionsBuffer.back()->att->significant_value == 1){
            ++run_1;
            this->order4_SubPartitionsBuffer.pop_back();
        }
        while(this->order4_SubPartitionsBuffer.size() > 0 && this->order4_SubPartitionsBuffer.back()->att->significant_value == 0) {
            ++run_0;
            this->order4_SubPartitionsBuffer.pop_back();
        }
        if(run_1 > 0){
            bits_run += this->BitsExpGolomb(run_1);
            run_1 = 0;
        }
        if(run_0 > 0){
            bits_run += this->BitsExpGolomb(run_0);
            run_0 = 0;
        }
    }
    this->file <<
       this->props->light_field_name << SEP <<
       this->props->hypercubo << SEP <<
       this->props->channel << SEP <<
       (8 + bits_run) << SEP << endl; // 8 bits do Last + bits referente a corrida
}

void Tree::ComputeMultipleLevesOfLastRun(){
    int run = 0;
    int bits_run = 0;
    int total_bits = 0;
    while (this->order8_SubPartitionsBuffer.size() > 0){
        if (this->order8_SubPartitionsBuffer.back()->att->significant_value != 0){
            total_bits += 4; // 4 bits para indicar uma posição em 16 subdivisões
            break;
        }else{
            this->order8_SubPartitionsBuffer.pop_back();
        }
    }
    for (int i = 0; i < this->order8_SubPartitionsBuffer.size(); ++i) { // adiciona os filhos das subpartisões restantes ao buffer de ordem 4
        for (int j = 0; j < this->order8_SubPartitionsBuffer[i]->child.size(); ++j) {
            this->order4_SubPartitionsBuffer.push_back(this->order8_SubPartitionsBuffer[i]->child[j]);
        }
    }
    this->order8_SubPartitionsBuffer.clear();
    while (this->order4_SubPartitionsBuffer.size() > 0){ // Encontra o Last no buffer de ordem 4
        if (this->order4_SubPartitionsBuffer.back()->att->significant_value != 0){
            total_bits += 8;
            break;
        }else{
            this->order4_SubPartitionsBuffer.pop_back();
        }
    }
    total_bits += this->order4_SubPartitionsBuffer.size();
    this->file <<
               this->props->light_field_name << SEP <<
               this->props->hypercubo << SEP <<
               this->props->channel << SEP <<
               (total_bits + bits_run) << SEP << endl; // bits dos Last + bits referente a corrida no nível ordem 4
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
/*    for (int i = 0; i < this->order4_SubPartitionsBuffer.size(); ++i) {
        cout << i << ": (x: " << this->order4_SubPartitionsBuffer[i]->node_pos.x <<
                ", y: " << this->order4_SubPartitionsBuffer[i]->node_pos.y <<
                ", u: " << this->order4_SubPartitionsBuffer[i]->node_pos.u <<
                ", v: " << this->order4_SubPartitionsBuffer[i]->node_pos.v << ")" << endl;
    }*/
}

int Tree::BitsExpGolomb(int code) {
    return (floor(log2(code + 1))) * 2 + 1;
}