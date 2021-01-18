#!/bin/bash

# Directories
#BINARY="./cmake-build-release/LF_Codec"
BINARY="./cmake-build-debug/LF_Codec"
DATASET_DIR="/home/igor/Documentos/Git/Full_datasets"

RESULT_DIR_ARITH="./results/Statistics/Arithmetic"
RESULT_DIR_LRE="./results/Statistics/LRE"

DATASETS=(Bikes DangerDeMort FountainVincent2 StonePillarsOutside)
LAMBDA=1
PREDICTION_MODE=none
TRANSFORM=(DCT_II)
QPS=(1 3 7 20)
NODES=0
ENTROPY_TYPE=(arithmetic lre)
LOG_OUTPUT=yes
FLAGS="-lytro -verbose -experimental"

function simulation() {
  JOINED_TRANSFORM=$(printf "__%s" "${TRANSFORM[@]}")
  JOINED_TRANSFORM=${JOINED_TRANSFORM:2}
  SIMULATION_ID="D${DATASET}_T${JOINED_TRANSFORM}_Q${QP}_X${QP}_Y${QP}_U${QP}_V${QP}_L${LAMBDA}_N${NODES}"

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
    -lfx 625 -lfy 434 -lfu 13 -lfv 13 \
    -blx 15 -bly 15 -blu 13 -blv 13 \
    -qx ${QP} -qy ${QP} -qu ${QP} -qv ${QP} -qp ${QP} \
    -lambda "${LAMBDA}" \
    -prediction "${PREDICTION_MODE}" \
    -entropy-type "${TYPE}" \
    -use-transforms ${TRANSFORM[@]} \
    -quadtree-max-inner-nodes ${NODES} \
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