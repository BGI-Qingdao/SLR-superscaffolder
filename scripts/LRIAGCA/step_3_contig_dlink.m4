#!/bin/bash

STEP="step_3 "

source ./__common_function.sh || exit 1

check_input xxx.cluster xxx.Arc xxx.updated.edge
try_backup_list xxx.connInfo
BIN/ContigDlink --prefix xxx --kvalue KVALUE --thread THREADS  --searchDepth SLEN 2>>log_contigdlink
check_output xxx.connInfo

