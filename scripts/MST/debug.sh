#!/bin/bash
source conf.ini
# prepare scripts
DATA="./"
MACRO=" -D xxx=$PREFIX \
   -D BIN=$EXEC_DIR \
   -D BWA=$BWA_DIR/bwa \
   -D MST_BIN_SIZE=$MST_BIN_SIZE \
   -D HT_BIN_SIZE=$HT_BIN_SIZE \
   -D BIN_THRESHOLD=$BIN_CLUSTER \
   -D HT_THRESHOLD=$HT_CLUSTER \
   -D MST_THRESHOLD=$MST_CLUSTER \
   -D MMCONTIG=$MIN_SCONTIG
   -D MIN_N=$MIN_N \
   -D MIN_C=$MIN_C \
   -D R1=$R1 \
   -D R2=$R2 \
   -D RANK=$RANK \
   -D SEED_MIN=$SEED_MIN \
   -D THREADS=$THREADS \
   -D BWA_K=$BWA_K"

echo "Generate all scripts ..."
m4 $MACRO $DATA/step_0_prepare_info.m4 >step_0_prepare_info.sh
m4 $MACRO $DATA/step_1_calc_seeds.m4 >step_1_calc_seeds.sh
m4 $MACRO $DATA/step_2_bin_cluster.m4 >step_2_bin_cluster.sh
m4 $MACRO $DATA/step_3_mst.m4 >step_3_mst.sh
m4 $MACRO $DATA/step_4_gap_oo.m4 >step_4_gap_oo.sh
m4 $MACRO $DATA/step_5_trunk2scaff.m4 >step_5_trunk2scaff.sh
m4 $MACRO $DATA/clean_prepare.m4 >clean_prepare.sh

chmod u+x *.sh

echo "prepare.sh finish"
