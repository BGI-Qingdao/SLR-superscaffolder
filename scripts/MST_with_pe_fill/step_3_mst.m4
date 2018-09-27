#!/bin/bash

STEP="step_3 "

source ./__common_function.sh || exit 1

check_input xxx.cluster
try_backup_list xxx.mintree xxx.mintree_trunk xxx.mintree_trunk_linear
BIN/MinTree --prefix xxx --threshold MST_THRESHOLD 2>>log_mst
check_output xxx.mintree xxx.mintree_trunk xxx.mintree_trunk_linear

