#!/bin/bash

STEP="step_5 "

source ./__common_function.sh || exit 1

check_input xxx.contig xxx.mintree_trunk_linear xxx.gap_oo 
try_backup_list xxx.scaff_seqs
BIN/Trunk2Scaff --prefix xxx  --K KVALUE  --min_scontig  MMCONTIG  2>>log_trunk2scaff
check_output xxx.scaff_seqs
