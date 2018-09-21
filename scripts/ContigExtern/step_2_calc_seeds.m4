#!/bin/bash

STEP="step_2 "
source ./__common_function.sh || exit 1


check_input xxx.contig 
try_backup_list xxx.seeds  xxx.pe_seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --min SEED_MIN  2>>log_staticsticunique
check_output xxx.seeds
mv xxx.seeds xxx.pe_seeds
check_output xxx.pe_seeds
