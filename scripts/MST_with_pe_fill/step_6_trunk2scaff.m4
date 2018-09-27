#!/bin/bash

STEP="step_6 "

source ./__common_function.sh || exit 1

check_input xxx.contig xxx.mintree_trunk_linear xxx.gap_oo xxx.trunk_fill
try_backup_list xxx.scaff_seqs
BIN/Trunk2Scaff --prefix xxx  --K KVALUE --gap_petrunk PE_FILL --gap_pe PE_FILL   --min_scontig  MMCONTIG  2>>log_trunk2scaff
check_output xxx.scaff_seqs
