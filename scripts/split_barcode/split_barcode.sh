#!/bin/bash

### 1st. Check parameters 

if [[ $# != 2 ]] ; then 
    echo "Usage : $0 raw_read1 raw_read2 "
    exit 1
fi

raw_read1=$1
raw_read2=$2
output_dir=$3

if [[ ! -f $raw_read1 || ! -f $raw_read2 ]] ; then 
    echo "Error : input raw reads not exist !!! exit ..."
    exit 1
fi

### 2nd. Check split_barcode_stLFR.pl && barcode_list.txt
base_dir=`dirname $0`
split_bc=$base_dir/split_barcode.pl
bc_list=$base_dir/barcode_list.txt
if [[ ! -f $split_bc || ! -f $bc_list ]] ; then 
    echo "Error : $split_bc or $bc_list is not exist !!! exit ... "
    exit 1
fi


### 3rd. Split barcodes
if [[ -f split_reads.1.fq.gz || -f split_reads.2.fq.gz ]] ; then 
    echo "split_reads.1.fq.gz or split_reads.2.fq.gz exist !!! exit ..."
    exit 1
fi
perl $split_bc $bc_list $raw_read1 $raw_read2 split_reads

