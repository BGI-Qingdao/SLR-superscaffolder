#!/usr/bin/env bash
if [[ $# != 3 ]] ; then
    echo "Usage: $0 contigOnRef barcodeOnRef outPrefix"
    exit
fi

contigOnRef=$1
barcodeOnRef=$2
outPrefix=$3
if [[ ! -f $contigOnRef ]] ; then
    echo "ERROR : $contigOnRef invalid ."
    exit
fi
if [[ ! -f $barcodeOnRef]] ; then
    echo "ERROR : $contigOnRef invalid ."
    exit
fi

barcodeOnContig=$outPrefix".barcodeOnRef"
barcodeOnBin=$outPrefix".barcodeOnBin"

./BarcodeOnContig -i $barcodeOnRef -c $contigOnRef -o $barcodeOnContig
./BarcodeOnBin -i $barcodeOnContig -o $barcodeOnBin -b 400




