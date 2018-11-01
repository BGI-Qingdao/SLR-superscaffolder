#!/bin/bash

STEP="step_5 "

source ./__common_function.sh || exit 1

check_input xxx.contigroad xxx.Arc xxx.updated.edge xxx.cluster xxx.barcodeOnContig
try_backup_list xxx.contigroadfill
BIN/FillContigRoad --prefix xxx  --kvalue KVALUE   --searchDepth SLEN  --fill_strategy 1  --thread THREADS --Ecov 15  2>>log_fillcontigroad
check_output xxx.contigroadfill
