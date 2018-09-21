#!/bin/bash

STEP="step_4 "

source ./__common_function.sh || exit 1


check_input chr19_soap2.contig chr19_soap2.Arc  chr19_soap2.updated.edge  chr19_soap2.contigroadfill
try_backup_list  chr19_soap2.super_used chr19_soap2.super_only \
                chr19_soap2.super_and_left chr19_soap2.contig_round_1\
                chr19_soap2.Arc_round_1 chr19_soap2.updated.edge_round_1 \
                chr19_soap2.ContigIndex_round_1

/home/guolidong/stLFR/debug/stLFR_Assembler/src/main/MergeContig --prefix chr19_soap2  --kvalue 63 --strategy 5 2>>log_mergecontig
check_output chr19_soap2.super_used chr19_soap2.super_only \
                chr19_soap2.super_and_left chr19_soap2.contig_round_1\
                chr19_soap2.Arc_round_1 chr19_soap2.updated.edge_round_1 \
                chr19_soap2.ContigIndex_round_1
