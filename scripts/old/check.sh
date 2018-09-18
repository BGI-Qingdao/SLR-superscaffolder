#!/bin/bash

# alian contig to ref
#bwa mem -k 63 -t 30 -S -p -a /research/rv-02/home/wangwenchao/ref/chr19/chr19.fa chr19_standard.contig >chr19_contig2ref.sam 
# get contigPos
#../bin/LinearClusterResult -i chr19_standard.cluster -p -c chr19_contig2ref.sam >chr19_contigPos.txt
# get linear seeds
#cat chr19_contigPos.txt | awk '{printf("%d\t%d\t%d\n", $4,$1,$2);}'|sort >sort_chr19_contigLinear.txt 
#cat chr19_seed.txt | sort >sort_chr19_seed.txt
#join sort_chr19_seed.txt sort_chr19_contigLinear.txt| awk '{printf("%d\t%d\t%d\t%d\n",$1,$2,$3,$4);}'  |  sort -nk 3 >chr19_seedLinear.txt
# align super to contigPos
../bin/LinkCheck -p chr19_seedLinear.txt <chr19_super.txt >check_log 2>&1 
