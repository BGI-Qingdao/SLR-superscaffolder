#!/bin/bash

STEP="step_5 "

source ./__common_function.sh || exit 1

check_input xxx.contig
try_backup_list xxx.gap.seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --middle_name gap --min GAP_BIN_SIZE 2>>log_staticsticunique
check_output xxx.gap.seeds

check_input xxx.gap.seeds xxx.barcodeOnContig
try_backup_list xxx.gap.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size GAP_BIN_SIZE --middle_name gap  --work_mode 2 --max_bin_size 3000000 2>>gap_chopbin
check_output xxx.gap.barcodeOnBin

check_input xxx.gap.barcodeOnBin  xxx.mintree_trunk_linear
try_backup_list xxx.gap.bin_cluster  xxx.gap.bin_cluster 
BIN/BinCluster --prefix xxx --thread THREADS  --middle_name gap  --threshold 0.01 --work_mode 1 --nb_only --bin_same_contig --pbc 2>>log_gapcluster
check_output xxx.gap.cluster xxx.gap.bin_cluster

check_input xxx.gap.bin_cluster  xxx.mintree_trunk_linear
try_backup_list xxx.gap_area xxx.gap_sim
BIN/GapSize --prefix xxx --max_gap 30000  2>>log_gap_size
check_output xxx.gap_area  xxx.gap_sim

