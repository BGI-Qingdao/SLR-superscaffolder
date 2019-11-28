#!/bin/bash

echo " Usage : ./MST_err_type.sh xxx.mst_error seed_types.txt"
ME=$1
echo " xxx.mst_error : $ME"
ST=$2
echo " seed_type.txt : $ST"

if [[ ! -e $ME || ! -e $ST ]] ; then 
    echo " invaid input , exit ..."
    exit 
fi
sort -k1b,1 $ME >$ME"_sorted"
sort -k1b,1 $ST >$ST"_sorted"
join $ME"_sorted" $ST"_sorted" | awk '{print $2}' |sort | uniq -c
