#!/bin/bash

echo " Usage : ./MST_linear_info_v1.sh xxx.mst_linear_info seed_types.txt xxx.ContigIndex"
LI=$1
echo " xxx.mst_linear_info : $LI"
ST=$2
echo " seed_type.txt : $ST"
CI=$3
echo " xxx.ContigIndex : $CI"
if [[ ! -e $LI || ! -e $ST || ! -e $CI ]] ; then 
    echo " invaid input , exit ..."
    exit 
fi
sort -k1b,1 $LI >$LI"_sorted"
sort -k1b,1 $ST >$ST"_sorted"
sort -k1b,1 $CI >$CI"_sorted"
echo ">mid type freq :"
join $LI"_sorted" $ST"_sorted" | awk '{print $4}' |sort | uniq -c

join $LI"_sorted" $ST"_sorted" >tmp_v1.txt
echo ""
join tmp_v1.txt $CI"_sorted" | awk 'print $2; }' | sort | uniq -c
echo ">non_unique common , f23:"
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "not_unique") printf("%s\t%s\n" $3,$4); }' >mst.non_unique.xy
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "not_unique") print $2; }' | sort | uniq -c

echo ">unique <=25K common barcode freq :"
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "is_unique" && ((0+$6) <= 25000) ) printf("%s\t%s\n" $3,$4); }' >mst.m25K.xy
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "is_unique" && ((0+$6) <= 25000) ) print $2; }' | sort | uniq -c

echo ">unique >25K common barcode freq :"
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "is_unique" && ((0+$6) > 25000) ) printf("%s\t%s\n" $3,$4); }' >mst.g25K.xy
join tmp_v1.txt $CI"_sorted" | awk '{if($5 == "is_unique" && ((0+$6) > 25000) ) print $2; }' | sort | uniq -c
