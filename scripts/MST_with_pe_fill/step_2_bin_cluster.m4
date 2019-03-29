#!/bin/bash

STEP="step_2 "

source ./__common_function.sh || exit 1

check_input xxx.seeds xxx.barcodeOnContig
try_backup_list xxx.mst.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size MST_BIN_SIZE  --work_mode 1 --max_bin_size 30000 --middle_name mst 2>>log_chopbin
check_output xxx.mst.barcodeOnBin

check_input xxx.mst.barcodeOnBin
try_backup_list xxx.mst.cluster
BIN/BinCluster --prefix xxx --thread THREADS --threshold MST_BIN_THRESHOLD --work_mode 1 --middle_name mst  2>>log_bincluster
check_output xxx.mst.cluster
