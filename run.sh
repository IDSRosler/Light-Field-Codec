#!/bin/bash

BINARY="./build/LF_Codec"
DATASET=~/lfcodec/datasets/Stone/
RESULT_DIR=../results-single-ish


function simulation {
    # gdb --args \
    # valgrind --leak-check=yes \
    DIR="$RESULT_DIR/$1.QP$QP/"
    mkdir -p $DIR
    $BINARY -input $DATASET -output $DIR -lytro \
        -blx 15 -bly 15 -blu 13 -blv 13 \
        -qx $QX -qy $QY -qu $QU -qv $QV -qp $QP \
        -lfx 625 -lfy 434 -lfu 13 -lfv 13 \
        -transform $1
}
QP=1
QX=1
QY=1
QU=1
QV=1
simulation MULTI 
