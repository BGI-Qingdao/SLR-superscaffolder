#!/bin/bash

mkdir output
mkdir logs
mkdir tmps

mv xxx.contig_round_1 output/xxx.contig || echo xxx.contig_round_1 is not exist!!!
mv xxx.Arc_round_1 output/xxx.Arc || echo xxx.Arc_round_1 is not exist!!! 
mv xxx.updated.edge_round_1 output/xxx.updated.edge || echo  xxx.updated.edge_round_1 is not exist!!! 
mv xxx.ContigIndex_round_1 output/xxx.ContigIndex || echo  xxx.ContigIndex_round_1 is not exist!!!·
mv xxx.super_only output/ || echo xxx.super_only is not exist!!!
mv xxx.super_and_left output/ || echo xxx.super_and_left is not exist!!!


mv log_* logs/
mv xxx.* tmps/
mv __start* tmps/
mv __end* tmps/

echo "info : output scaff is in output/scaff_seqs "
echo "info : logs are all in logs/*"
echo "info : you can delete folder tmps by : rm -rf tmps "
