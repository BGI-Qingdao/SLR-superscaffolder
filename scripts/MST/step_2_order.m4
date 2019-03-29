#!/bin/bash

STEP="step_2 "

source ./__common_function.sh || exit 1

check_input xxx.contig
try_backup_list xxx.mst.seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --middle_name mst  --min MST_SEED_MIN 2>>log_staticsticunique
check_output xxx.mst.seeds

check_input xxx.mst.seeds xxx.barcodeOnContig
try_backup_list xxx.mst.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size MST_BIN_SIZE --middle_name mst  --work_mode 1 --max_bin_size 30000 2>>log_chopbin
check_output xxx.mst.barcodeOnBin

check_input xxx.mst.barcodeOnBin
try_backup_list xxx.mst.cluster
BIN/BinCluster --prefix xxx --thread THREADS  --middle_name mst  --threshold MST_BIN_THRESHOLD --work_mode 1 2>>log_bincluster
check_output xxx.mst.cluster

check_input xxx.cluster
try_backup_list xxx.mintree_trunk_linear
BIN/MST --prefix xxx --threshold MST_THRESHOLD --del_fac 0.95 --del_round 100 2>>log_mst
check_output xxx.mintree_trunk_linear

