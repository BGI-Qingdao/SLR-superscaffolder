#!/bin/bash

source ./__common_function.sh || exit 1

check_input xxx.seeds xxx.barcodeOnContig
try_backup_list xxx.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size BIN_SIZE  --delete_tail TAIL  2>>log_chopbin
check_output xxx.barcodeOnBin

check_input xxx.barcodeOnBin
try_backup_list xxx.cluster
BIN/BinCluster --prefix xxx --thread THREADS --threshold BIN_THRESHOLD 2>>log_bincluster
check_output xxx.cluster
