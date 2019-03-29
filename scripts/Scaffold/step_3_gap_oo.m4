#!/bin/bash

STEP="step_3 "

source ./__common_function.sh || exit 1

check_input xxx.contig
try_backup_list xxx.oo.seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --middle_name  oo  --min MST_BIN_SIZE 2>>log_staticsticunique
check_output xxx.oo.seeds

check_input xxx.oo.seeds xxx.barcodeOnContig
try_backup_list xxx.oo.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size HT_BIN_SIZE --middle_name oo  --work_mode 2  2>>log_chopbin_ht
check_output xxx.oo.barcodeOnBin

check_input xxx.oo.barcodeOnBin
try_backup_list xxx.oo.cluster xxx.oo.bin_cluster
BIN/BinCluster --prefix xxx --thread THREADS --middle_name oo  --work_mode 1 --threshold HT_BIN_THRESHOLD --pbc --bin_same_contig  2>>log_bincluster_ht
check_output xxx.oo.cluster xxx.oo.bin_cluster

check_input xxx.mintree_trunk_linear xxx.oo.bin_cluster
try_backup_list xxx.gap_oo
BIN/Orientation --prefix xxx --rank RANK   2>>log_gap_oo
check_output xxx.gap_oo
