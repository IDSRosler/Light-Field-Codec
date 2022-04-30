#!/bin/bash

# Directories
BINARY="./cmake-build-release/LF_Codec"
#BINARY="./cmake-build-debug/LF_Codec"
DATASET_DIR="/media/idsrosler/67b63062-c831-4d62-8fbe-11385746955f/igor/Git/CTC_Datasets/Lenslets/"

#RESULT_DIR_ARITH="/media/idsrosler/67b63062-c831-4d62-8fbe-11385746955f/igor/Git/TCC-output/Arithmetic/Lenslet/"
RESULT_DIR_LRE="/media/idsrosler/67b63062-c831-4d62-8fbe-11385746955f/igor/Git/TCC-output/RLE/Lenslet/"

DATASETS=(Bikes Danger_de_Mort Fountain_Vincent2 Stone_Pillars_Outside)
LAMBDA=1
PREDICTION_MODE=none
TRANSFORM=(DCT_II)
SPATIAL_SIZE=15
ANGULAR_SIZE=15
QPS=(1 3 7 20)
QPX=1
QPY=1
QPU=1
QPV=1
NODES=0
ENTROPY_TYPE=(lre)
LOG_OUTPUT=yes
FLAGS="-lytro -verbose -experimental -export-statistics"

function simulation() {
  JOINED_TRANSFORM=$(printf "__%s" "${TRANSFORM[@]}")
  JOINED_TRANSFORM=${JOINED_TRANSFORM:2}
  SIMULATION_ID="D${DATASET}_T${JOINED_TRANSFORM}_Q${QP}_X${QPX}_Y${QPY}_U${QPU}_V${QPV}_L${LAMBDA}_N${NODES}"

  mkdir -p "${RESULT}/${SIMULATION_ID}/"
  mkdir -p "${RESULT}/logs/"

  echo "Running ${SIMULATION_ID}"
  if [ $LOG_OUTPUT == "yes" ]; then
    LOGFILE="${RESULT}/logs/${SIMULATION_ID}.log"
    FLAGS="${FLAGS} -verbose"
  else
    LOGFILE="/dev/null"
    FLAGS="${FLAGS} -show-progress-bar"
  fi

  $BINARY \
    -input "${DATASET_DIR}/${DATASET}/" \
    -output "${RESULT}/${SIMULATION_ID}/" \
    -lfx 625 -lfy 434 -lfu 15 -lfv 15 \
    -blx 15 -bly 15 -blu 15 -blv 15 \
    -qx ${QPX} -qy ${QPY} -qu ${QPU} -qv ${QPV} -qp ${QP} \
    -lambda "${LAMBDA}" \
    -prediction "${PREDICTION_MODE}" \
    -entropy-type "${TYPE}" \
    -use-transforms ${TRANSFORM[@]} \
    -partition-tree-max-depth ${NODES} \
    -transform-min-spatial-size ${SPATIAL_SIZE} \
    -transform-min-angular-size ${ANGULAR_SIZE} \
    ${FLAGS} | tee "${LOGFILE}"
}

for TYPE in ${ENTROPY_TYPE[@]};
do
  for DATASET in ${DATASETS[@]};
  do
    for QP in ${QPS[@]};
    do
      if [ "$TYPE" == "arithmetic" ];
      then
        RESULT=$RESULT_DIR_ARITH
      else
        RESULT=$RESULT_DIR_LRE
      fi
      simulation
    done
  done
done