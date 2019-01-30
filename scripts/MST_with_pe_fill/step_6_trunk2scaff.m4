#!/bin/bash
STEP="step_6 "

source ./__common_function.sh || exit 1

check_input xxx.contig xxx.mintree_trunk_linear xxx.gap_oo xxx.gap_area xxx.trunk_fill
try_backup_list xxx.scaff_infos
BIN/Trunk2ScaffInfo --prefix xxx 2>>log_trunk2scaffinfo
check_output xxx.scaff_infos

check_input xxx.contig xxx.scaff_infos
try_backup_list xxx.scaff_seqs
BIN/ScaffInfo2Seq --prefix xxx --min_n MIN_N --min_c MIN_C 2>>log_scaffinfo2seq
check_output xxx.scaff_seqs
