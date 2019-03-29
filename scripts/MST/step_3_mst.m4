#!/bin/bash

STEP="step_3 "

source ./__common_function.sh || exit 1

check_input xxx.cluster
try_backup_list xxx.mintree_trunk_linear
BIN/MST --prefix xxx --threshold MST_THRESHOLD --del_fac 0.95 --del_round 100 2>>log_mst
check_output xxx.mintree_trunk_linear

