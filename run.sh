#!/bin/bash

# Default values
BINARY="./cmake-build-release/LF_Codec"
DATASET_DIR="../LFCodecResources/Datasets"
RESULT_DIR="./results"

QP=1
QX=1
QY=1
QU=1
QV=1
LAMBDA=1
DATASET=Bikes
TRANSFORM=DCT_II
LOG_OUTPUT=no
FLAGS="-lytro -experimental"

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
    -use-transforms ${TRANSFORM[@]} \
    -quadtree-max-inner-nodes ${NODES} \
    ${FLAGS} | tee "${LOGFILE}"
}

LOG_OUTPUT=yes

TRANSFORM=(DCT_II)
NODES=4
FLAGS="-lytro -verbose -experimental"
simulation
#
#FLAGS="${FLAGS} -experimental"
#simulation
