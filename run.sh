#!/bin/bash

# Default values
if [ "$BINARY" == "" ]; then
  echo "Variable BINARY not found. Using default value."
  BINARY="./cmake-build-debug/LF_Codec"
  #BINARY="./cmake-build-release/LF_Codec"
fi
if [ "$DATASET_DIR" == "" ]; then
  echo "Variable DATASET_DIR not found. Using default value."
  DATASET_DIR="/home/italo/Documentos/Light-Field-Codec/LFs/"
fi
if [ "$RESULT_DIR" == "" ]; then
  echo "Variable RESULT_DIR not found. Using default value."
  RESULT_DIR="./results/test/LRE"
fi

QP=20
QX=20
QY=20
QU=20
QV=20
LAMBDA=1
DATASET=Bikes
TRANSFORM=DCT_II
PREDICTION_MODE=none
ENTROPY_TYPE=lre
LOG_OUTPUT=no
FLAGS="-lytro -experimental -lossless"

function simulation() {
  JOINED_TRANSFORM=$(printf "__%s" "${TRANSFORM[@]}")
  JOINED_TRANSFORM=${JOINED_TRANSFORM:2}
  SIMULATION_ID="D${DATASET}_T${JOINED_TRANSFORM}_Q${QP}_X${QX}_Y${QY}_U${QU}_V${QV}_L${LAMBDA}_N${NODES}"

  mkdir -p "${RESULT_DIR}/${SIMULATION_ID}/"
  mkdir -p "${RESULT_DIR}/logs/"

  echo "Running ${SIMULATION_ID}"
  if [ $LOG_OUTPUT == "yes" ]; then
    LOGFILE="${RESULT_DIR}/logs/${SIMULATION_ID}.log"
    FLAGS="${FLAGS} -verbose"
  else
    LOGFILE="/dev/null"
    FLAGS="${FLAGS} -show-progress-bar"
  fi

  $BINARY \
    -input "${DATASET_DIR}/${DATASET}/" \
    -output "${RESULT_DIR}/${SIMULATION_ID}/" \
    -lfx 625 -lfy 434 -lfu 13 -lfv 13 \
    -blx 15 -bly 15 -blu 13 -blv 13 \
    -qx ${QX} -qy ${QY} -qu ${QU} -qv ${QV} -qp ${QP} \
    -lambda "${LAMBDA}" \
    -prediction "${PREDICTION_MODE}" \
    -entropy-type "${ENTROPY_TYPE}" \
    -use-transforms ${TRANSFORM[@]} \
    -quadtree-max-inner-nodes ${NODES} \
    ${FLAGS} | tee "${LOGFILE}"
}

LOG_OUTPUT=yes

TRANSFORM=(DCT_II)
NODES=0
FLAGS="-lytro -verbose -export-statistics -"
simulation
#
#FLAGS="${FLAGS} -experimental"
#simulation
