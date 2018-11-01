#!/bin/bash

STEP="step_3 "
source ./__common_function.sh || exit 1

check_input xxx.pe_seeds  xxx.updated.edge  xxx.Arc xxx.barcode_at_contig_v1 xxx.contig_pe_conns
try_backup_list xxx.seed_extern_fill xxx.contigroadfill
BIN/ExternContigByPE  --prefix xxx --kvalue KVALUE  --min_count PE_MIN --min_factor PE_FAC --min_bcount BC_MIN --min_bfactor BC_FAC  2>>log_externcontig
check_output  xxx.seed_extern_fill
cp  xxx.seed_extern_fill  xxx.contigroadfill
check_output xxx.contigroadfill

