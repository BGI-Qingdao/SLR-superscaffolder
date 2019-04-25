#!/bin/bash
#source conf.ini
# check parameters
if [[ $# != 1 ]] ; then
    echo "Usage : $0 conf_file"
    exit 1
fi

# check parameters as file
if [[  -f "$1" ]] ; then 
    source "$1"
else
    echo "file $1 is not exsist !!! exit ... "
    exit 1
fi

# check install dir
INSTALL_DIR=$STLFR_ASSEMBLER_DIR
EXEC_DIR=$INSTALL_DIR/bin
DATA=$INSTALL_DIR/scaffold

# check target directory
if [[ ! -d $EXEC_DIR || ! -e $EXEC_DIR/FakeSOAPContig ]] ; then 
    echo $EXEC_DIR is not a valid stLFR Assembler bin dir . exit ..."
    exit 1
fi
if [[ ! -d $DATA || ! -f $DATA/step_1_prepare_info.m4 ]] ; then 
    echo $DATA is not a valid stLFR Assembler Scaffold dir . exit ..."
    exit 1
fi

# check SOAP diectory
SOAP=$SOAP_DIR
Contig=$SOAP/$PREFIX".contig"
ContigIndex=$SOAP/$PREFIX".ContigIndex"
if [[ ! -e $Contig || ! -e $ContigIndex  ]] ; then 
    echo "$SOAP/$PREFIX is not a valid SOAP prefix . exit ..."
    exit 1
fi

# check target directory
if [[ -d $PROJECT_NAME ]] ; then 
    echo "Project $PROJECT_NAME already exsist !!! exit ..."
    echo "you can enter directory $PROJECT_NAME and run commands directly ."
    exit 1
fi

# mkdir target file
echo "mkdir $PROJECT_NAME"
mkdir $PROJECT_NAME
# ln input file
cd $PROJECT_NAME
echo "ln -s $Contig"
ln -s $Contig
echo "ln -s $ContigIndex"
ln -s $ContigIndex
chmod a-w $PREFIX*
chmod a-x $PREFIX*
# prepare scripts

MACRO=" -D xxx=$PREFIX \
   -D BWA=$BWA_DIR/bwa \
   -D BIN=$EXEC_DIR \
   -D HT_BIN_SIZE=$HT_BIN_SIZE \
   -D HT_BIN_THRESHOLD=$HT_BIN_CLUSTER \
   -D KVALUE=$SOAP_K \
   -D MST_BIN_SIZE=$MST_BIN_SIZE \
   -D MST_BIN_THRESHOLD=$MST_BIN_CLUSTER \
   -D MST_THRESHOLD=$MST_CLUSTER \
   -D MMCONTIG=$MIN_SCONTIG \
   -D MIN_N=$MIN_N \
   -D GAP_BIN_SIZE=$GAP_BIN_SIZE \
   -D MIN_C=$MIN_C \
   -D R1=$R1 \
   -D R2=$R2 \
   -D RANK=$RANK \
   -D THREADS=$THREADS \
   -D MAX_IS=$MAX_INSERT_SIZE\
   -D PE_SEARCH_MAX=$PE_SEARCH_MAX \
   -D PE_SEED_MIN=$PE_SEED_MIN \
   -D PE_MIN_B=$PE_MIN_JOINBARCODES \
   -D PE_MIN_COUNT=$PE_MIN_COUNT \
   -D PE_FILL=$PE_FILL \
   -D BWA_K=$BWA_K"

echo "Generate all scripts ..."
cp $DATA/__common_function.sh ./
cp $DATA/run.sh run.sh
m4 $MACRO $DATA/step_1_prepare_info.m4 >step_1_prepare_info.sh
m4 $MACRO $DATA/step_2_order.m4 >step_2_order.sh
m4 $MACRO $DATA/step_3_orientation.m4 >step_3_orientation.sh
m4 $MACRO $DATA/step_4_pe_fill.m4 >step_4_pe_fill.sh
m4 $MACRO $DATA/step_5_gapsize.m4 >step_5_gapsize.sh
m4 $MACRO $DATA/step_6_gen_seq.m4 >step_6_gen_seq.sh
m4 $MACRO $DATA/clean_prepare.m4 >clean_prepare.sh

chmod u+x *.sh 

echo "prepare.sh finish"
