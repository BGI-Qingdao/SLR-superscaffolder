#!/bin/bash

mkdir output
mkdir logs
mkdir tmps

mv xxx.scaff_infos output/ || echo xxx.scaff_infos is not exist!!!
mv xxx.scaff_seqs output/ || echo xxx.scaff_seqs is not exist!!!
mv xxx.scaff_agp output/ || echo xxx.scaff_seqs is not exist!!!

mv log_* logs/
mv xxx.* tmps/
mv __start* tmps/
mv __end* tmps/

echo "info : output scaff is in output/"
echo "info : logs are all in logs/*"
echo "info : you can delete folder tmps by : rm -rf tmps "
