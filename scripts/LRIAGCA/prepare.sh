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

# check SOAP diectory
SOAP=$SOAP_DIR
Contig=$SOAP/$PREFIX".contig"
Edge=$SOAP/$PREFIX".updated.edge"
Arc=$SOAP/$PREFIX".Arc"
Index=$SOAP/$PREFIX".ContigIndex"
if [[ ! -e $Contig ||  ! -e $Edge || ! -e $Arc || ! -e $Index ]] ; then 
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
echo "ln -s $Edge"
ln -s $Edge
echo "ln -s $Index"
ln -s $Index
echo "ln -s $Arc"
ln -s $Arc
chmod a-w $PREFIX*
chmod a-x $PREFIX*
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
cp $DATA/__common_function.sh ./
m4 $MACRO $DATA/step_0_prepare_info.m4 >step_0_prepare_info.sh
m4 $MACRO $DATA/step_1_calc_seeds.m4 >step_1_calc_seeds.sh
m4 $MACRO $DATA/step_2_bin_cluster.m4 >step_2_bin_cluster.sh
m4 $MACRO $DATA/step_3_contig_dlink.m4 >step_3_contig_dlink.sh
m4 $MACRO $DATA/step_4_contig_road.m4 >step_4_contig_road.sh
m4 $MACRO $DATA/step_5_contig_road_fill.m4 >step_5_contig_road_fill.sh
m4 $MACRO $DATA/step_6_merge_contig.m4 >step_6_merge_contig.sh
m4 $MACRO $DATA/run.m4 >run.sh
chmod u+x step_*.sh 
chmod u+x run.sh

echo "prepare.sh done"
