#!/bin/bash
source conf.ini
# prepare scripts

MACRO=" -D xxx=$PREFIX \
   -D BIN=$EXEC_DIR \
   -D BIN_SIZE=$BIN_SIZE \
   -D BWA=$BWA_DIR/bwa \
   -D KVALUE=$SOAP_K \
   -D R1=$R1 \
   -D R2=$R2 \
   -D SLEN=$SEARCH_DEPTH \
   -D SF=$SIM_FACTOR \
   -D THRESHOLD=$CLUSTER \
   -D THREADS=$THREADS \
   -D TAIL=$TAIL_DEL \
   -D LF=$LEN_FACTOR "

echo "Generate all scripts ..."
m4 $MACRO $DATA/step_0_prepare_info.m4 >step_0_prepare_info.sh
m4 $MACRO $DATA/step_1_calc_seeds.m4 >step_1_calc_seeds.sh
m4 $MACRO $DATA/step_2_bin_cluster.m4 >step_2_bin_cluster.sh
m4 $MACRO $DATA/step_3_contig_dlink.m4 >step_3_contig_dlink.sh
m4 $MACRO $DATA/step_4_contig_road.m4 >step_4_contig_road.sh
m4 $MACRO $DATA/step_5_contig_road_fill.m4 >step_5_contig_road_fill.sh
m4 $MACRO $DATA/step_6_merge_contig.m4 >step_6_merge_contig.sh
m4 $MACRO $DATA/run.m4 >run.sh
m4 $MACRO $DATA/clean_prepare.m4 >clean_prepare.sh

chmod u+x step_*.sh 
chmod u+x run.sh

echo "prepare.sh done"
