#!/bin/bash

### 1st. Check parameters 

if [[ $# != 3 ]] ; then 
    echo "Usage : $0 raw_read1 raw_read2 output_dir"
    exit 1
fi

raw_read1=$1
raw_read2=$2
output_dir=$3

if [[ ! -f $raw_read1 || ! -f $raw_read2 ]] ; then 
    echo "Error : input raw reads not exist !!! exit ..."
    exit 1
fi

if [[ ! -d $output_dir ]] ; then 
    echo "Info : $output_dir is not exist , create it now ..."
    mkdir $output_dir 
    if [[ ! -d $output_dir ]] ; then 
        echo "Error : failed to create $output_dir !!! exit ..."
        exit 1 ;
    fi
fi

### 2nd. Check split_barcode_stLFR.pl && barcode_list.txt
base_dir=`dirname $0`
split_bc=$base_dir/split_barcode_stLFR.pl
bc_list=$base_dir/barcode_list.txt
if [[ ! -f $split_bc || ! -f $bc_list ]] ; then 
    echo "Error : $split_bc or $bc_list is not exist !!! exit ... "
    exit 1
fi

### 3rd. Check sequence length of read2
len_r2=`less $raw_read2 | head | sed -n '2p' | awk '{print  length($1);}'`

if [[ $len_r2 == 142 ]] ; then 
    $split_bc -i1 $raw_read1 -i2 $raw_read2 -o $output_dir -b $bc_list -r "101_10,117_10,133_10"
elif [[ $len_r2 == 154 ]] ; then 
    $split_bc -i1 $raw_read1 -i2 $raw_read2 -o $output_dir -b $bc_list -r "101_10,117_10,145_10"
else
    echo "Error : unknow read2 sequence length : $len_r2 !!! exit ..."
    exit 1
fi
