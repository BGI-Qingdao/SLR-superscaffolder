#!/bin/bash

STEP="step_4 "

source ./__common_function.sh || exit 1

check_input xxx.seeds xxx.barcodeOnContig
try_backup_list xxx.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size HT_BIN_SIZE --ht_only  2>>log_chopbin_ht
check_output xxx.barcodeOnBin

check_input xxx.barcodeOnBin
try_backup_list xxx.cluster xxx.bin_cluster
BIN/BinCluster --prefix xxx --thread THREADS --threshold HT_THRESHOLD --pbc --bin_same_contig  2>>log_bincluster_ht
check_output xxx.cluster xxx.bin_cluster

check_input xxx.mintree_trunk_linear xxx.bin_cluster
try_backup_list xxx.gap_oo xxx.gap_area
BIN/Gap_OAO --prefix xxx --rank RANK  --calc_linear 2>>log_gap_oo
check_output xxx.gap_oo xxx.gap_area
