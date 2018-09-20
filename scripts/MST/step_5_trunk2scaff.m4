#!/bin/bash

source ./__common_function.sh || exit 1

check_input xxx.contig xxx.mintree_trunk_linear xxx.gap_oo 
try_backup_list xxx.scaff_seqs
BIN/FillContigRoad --prefix xxx  --kvalue KVALUE   --searchDepth SLEN  --fill_strategy 1  --thread THREADS --Ecov 15  2>>log_trunk2scaff
check_output xxx.scaff_seqs
