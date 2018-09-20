#!/bin/bash

STEP="step_1 "

source ./__common_function.sh || exit 1

check_input xxx.contig
try_backup_list xxx.seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --min eval(BIN_SIZE+TAIL) 2>>log_staticsticunique
check_output xxx.seeds
