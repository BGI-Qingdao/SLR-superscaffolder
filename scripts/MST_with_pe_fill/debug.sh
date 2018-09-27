#!/bin/bash
source conf.ini
# prepare scripts
DATA="./"

MACRO=" -D xxx=$PREFIX \
   -D BWA=$BWA_DIR/bwa \
   -D BIN=$EXEC_DIR \
   -D MST_BIN_SIZE=$MST_BIN_SIZE \
   -D HT_BIN_SIZE=$HT_BIN_SIZE \
   -D PE_BIN_SIZE=$PE_BIN_SIZE \
   -D MST_BIN_THRESHOLD=$MST_BIN_CLUSTER \
   -D MAX_IS=$MAX_INSERT_SIZE\
   -D HT_BIN_THRESHOLD=$HT_BIN_CLUSTER \
   -D MST_THRESHOLD=$MST_CLUSTER \
   -D MMCONTIG=$MIN_SCONTIG
   -D KVALUE=$SOAP_K \
   -D PE_BIN_THRESHOLD=$PE_BIN_CLUSTER \
   -D PE_SEARCH_MAX=$PE_SEARCH_MAX \
   -D PE_MIN_COUNT=$PE_MIN_COUNT \
   -D PE_FILL=$PE_FILL \
   -D R1=$R1 \
   -D R2=$R2 \
   -D RANK=$RANK \
   -D MST_SEED_MIN=$MST_SEED_MIN \
   -D PE_SEED_MIN=$PE_SEED_MIN \
   -D THREADS=$THREADS \
   -D TAIL=$TAIL_DEL "

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
