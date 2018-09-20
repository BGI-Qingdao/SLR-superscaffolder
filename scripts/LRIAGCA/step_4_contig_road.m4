#!/bin/bash

STEP="step_4 "

source ./__common_function.sh || exit 1

check_input xxx.connInfo
try_backup_list xxx.contigroad
BIN/LinearCDG --prefix xxx --len_factor LF --sim_factor SF 2>>log_linearcdg
check_output xxx.contigroad
