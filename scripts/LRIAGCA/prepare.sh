#!/bin/bash
source conf.ini
## check parameters
#if [[ $# != 1 ]] ; then
#    echo "Usage : $0 conf_file"
#    exit 1
#fi
#
## check parameters as file
#if [[  -f "$1" ]] ; then 
#    source "$1"
#else
#    echo "file $1 is not exsist !!! exit ... "
#    exit 1
#fi
#
## check SOAP diectory
#Contig=$SOAP/$PREFIX".contig"
#Edge=$SOAP/$PREFIX".updated.edge"
#Arc=$SOAP/$PREFIX".Arc"
#Index=$SOAP/$PREFIX".ContigIndex"
#if [[ ! -e $Contig ||  ! -e $Edge || ! -e $Arc || ! -e $Index ]] ; then 
#    echo "$SOAP/$PREFIX is not a valid SOAP prefix . exit ..."
#    exit 1
#fi
#
## check target directory
#if [[ -d $PROJECT_NAME ]] ; then 
#    echo "Project $PROJECT_NAME already exsist !!! exit ..."
#    echo "you can enter directory $PROJECT_NAME and run commands directly ."
#    exit 1
#fi
#
## mkdir target file
#echo "mkdir $PROJECT_NAME"
#mkdir $PROJECT_NAME
## ln input file
#cd $PROJECT_NAME
#echo "ln -s $Contig"
#ln -s $Contig
#echo "ln -s $Edge"
#ln -s $Edge
#echo "ln -s $Index"
#ln -s $Index
#echo "ln -s $Arc"
#ln -s $Arc
#
# prepare scripts
m4 -D xxx=$PREFIX \
   -D THREADS=$THREADS\
   -D R1=$R1\
   -D R2=$R2\
   -D BIN=$BIN_DIR \
   $DATA/step_0_prepare_info.m4 >step_0_prepare_info.sh

