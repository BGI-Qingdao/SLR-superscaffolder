#!/bin/bash

STEP="step_2 "

source ./__common_function.sh || exit 1

check_input xxx.seeds xxx.barcodeOnContig
try_backup_list xxx.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size MST_BIN_SIZE --middle_name mst  --work_mode 1 --max_bin_size 30000 2>>log_chopbin
check_output xxx.barcodeOnBin

check_input xxx.barcodeOnBin
try_backup_list xxx.cluster
BIN/BinCluster --prefix xxx --thread THREADS  --middle_name mst  --threshold MST_BIN_THRESHOLD --work_mode 1 2>>log_bincluster
check_output xxx.cluster
