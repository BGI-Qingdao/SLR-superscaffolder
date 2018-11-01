#!/bin/bash

STEP="step_1 "

source ./__common_function.sh || exit 1

check_input xxx.contig 
try_backup_list xxx.contig_long xxx.contig_short
BIN/SplitContig --prefix xxx   --threshold READ_LEN 2>>log_split_read || exit 1
check_output xxx.contig_long xxx.contig_short

BWA index xxx.contig_long >log_bwa_index_long 2>&1 || exit 1
BWA index R1 >log_bwa_index_r1 2>&1 || exit 1
BWA index R2 >log_bwa_index_r2 2>&1 || exit 1

check_input xxx.contig_long R1 R2
try_backup_list xxx.read2contig.sam 
BWA mem -t THREADS -k eval(TAIL+1) xxx.contig_long R1 R2 >xxx.read2contig.sam 2>>log_bwa_mem_long ||exit 1
check_output xxx.read2contig.sam

check_input xxx.contig_short R1
try_backup_list xxx.contig2r1.sam
BWA mem -t THREADS  -k eval(TAIL+1)  -S -p -a  R1  xxx.contig_short >xxx.contig2r1.sam 2>>log_bwa_memr1 ||exit 1
check_output xxx.read2contig.sam 

check_input xxx.contig_short R2
try_backup_list xxx.contig2r2.sam 
BWA mem -t THREADS  -k eval(TAIL+1)  -S -p -a  R2  xxx.contig_short >xxx.contig2r2.sam 2>>log_bwa_memr2 ||exit 1
check_output xxx.contig2r2.sam 

check_input xxx.contig2r1.sam xxx.contig2r2.sam xxx.readNameList xxx.barcodeList
try_backup_list xxx.contig2read_v1
BIN/ParseContig2Read --prefix xxx 2>>log_parse_contig2read || exit 1
check_output xxx.contig2read_v1

check_input xxx.read2contig.sam  xxx.readNameList xxx.barcodeList
try_backup_list xxx.read2contig_v1
BIN/ParseRead2Contig --prefix xxx 2>>log_parse_read2contig || exit 1
check_output xxx.read2contig_v1

check_input xxx.read2contig_v1 xxx.contig2read_v1
try_backup_list xxx.contig_pe_conns xxx.barcode_at_contig_v1 xxx.contig_at_barcode_v1
BIN/MergePEInfo --prefix xxx 2>>log_mergepeinfo || exit 1 
check_output xxx.contig_pe_conns xxx.barcode_at_contig_v1 xxx.contig_at_barcode_v1

