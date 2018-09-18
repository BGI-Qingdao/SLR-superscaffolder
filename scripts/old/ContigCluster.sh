#!/bin/bash

prefix=""
threshold="0.03"
n=0
while getopts "o:s:n:" arg
do
    case $arg in 
        o)
            prefix=$optarg
            ;;
        s)
            threshold=$optarg
            ;;
        n)
            n=$optarg
            ;;
        ?)
            echo "Warning unknow $arg"
    esac
done

barcodeOnBin=$prefix".barcodeOnBin"
if [[ ! -f  $barcodeOnBin ]] ; then
    echo "$barcodeOnBin not exsist!"
    exit;
fi
if [[ n -lt 1 ]] ; then 
    n=`wc -l $barcodeOnBin`
fi

cluster_split=$prefix".cluster_split"
./BinCluster -i $barcodeOnBin -o  $cluster_split -s $threshold -n $n >"temp.$cluster_split"

