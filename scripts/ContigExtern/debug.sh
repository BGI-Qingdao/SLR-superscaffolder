#!/bin/bash
source conf.ini
# prepare scripts

DATA="./"

MACRO=" -D xxx=$PREFIX \
   -D BWA=$BWA_DIR/bwa \
   -D BC_MIN=$BC_MIN \
   -D BC_FAC=$BC_FAC \
   -D BIN=$EXEC_DIR \
   -D SEED_MIN=$SEED_MIN
   -D KVALUE=$SOAP_K \
   -D R1=$R1 \
   -D R2=$R2 \
   -D READ_LEN=$READ_LEN \
   -D TAIL=$TAIL_DEL \
   -D THREADS=$THREADS \
   -D PE_MIN=$PE_MIN \
   -D PE_FAC=$PE_FAC "

echo "Generate all scripts ..."
m4 $MACRO $DATA/step_0_parse_read1.m4 >step_0_parse_read1.sh
m4 $MACRO $DATA/step_1_prepare_info.m4 >step_1_prepare_info.sh
m4 $MACRO $DATA/step_2_calc_seeds.m4 >step_2_calc_seeds.sh
m4 $MACRO $DATA/step_3_extern_contig.m4 >step_3_extern_contig.sh
m4 $MACRO $DATA/step_4_merge_contig.m4 >step_4_merge_contig.sh
m4 $MACRO $DATA/clean_prepare.m4 >clean_prepare.sh

chmod u+x *.sh 
