#!/bin/bash

STEP="step_4 "

source ./__common_function.sh || exit 1

check_input xxx.contig xxx.Arc  xxx.updated.edge  xxx.contigroadfill
try_backup_list  xxx.super_used xxx.super_only \
                xxx.super_and_left xxx.contig_round_1\
                xxx.Arc_round_1 xxx.updated.edge_round_1 \
                xxx.ContigIndex_round_1

BIN/MergeContig --prefix xxx  --kvalue KVALUE --strategy 5 2>>log_mergecontig
check_output xxx.super_used xxx.super_only \
                xxx.super_and_left xxx.contig_round_1\
                xxx.Arc_round_1 xxx.updated.edge_round_1 \
                xxx.ContigIndex_round_1
