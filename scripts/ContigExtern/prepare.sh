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
ContigIndex=$SOAP/$PREFIX".ContigIndex"
if [[ ! -e $Contig || ! -e $ContigIndex ]] ; then 
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
cp $DATA/__common_function.sh ./
cp $DATA/run.sh ./
m4 $MACRO $DATA/step_0_parse_read1.m4 >step_0_parse_read1.sh
m4 $MACRO $DATA/step_1_prepare_info.m4 >step_1_prepare_info.sh
m4 $MACRO $DATA/step_2_calc_seeds.m4 >step_2_calc_seeds.sh
m4 $MACRO $DATA/step_3_extern_contig.m4 >step_3_extern_contig.sh
m4 $MACRO $DATA/step_4_merge_contig.m4 >step_4_merge_contig.sh
m4 $MACRO $DATA/clean_prepare.m4 >clean_prepare.sh

chmod u+x *.sh 

echo "prepare.sh finish"
