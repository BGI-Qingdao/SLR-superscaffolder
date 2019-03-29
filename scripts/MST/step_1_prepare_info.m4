#!/bin/bash

STEP="step_1 "

source ./__common_function.sh || exit 1

check_input xxx.contig R1 R2

BWA index xxx.contig >log_bwa_index 2>&1 || exit 1

try_backup_list xxx.read2contig.sam log_bwa_mem
check_input xxx.contig R1 R2
BWA mem -t THREADS -k BWA_K xxx.contig R1 R2 >xxx.read2contig.sam 2>>log_bwa_mem ||exit 1
check_output xxx.read2contig.sam log_bwa_mem

try_backup_list   xxx.barcodeList   xxx.readNameList
check_input R1
BIN/ParseReadName --prefix xxx --read1 R1 2>>log_parsereadname || exit 1
check_output xxx.barcodeList xxx.readNameList

try_backup_list xxx.read2contig
check_input xxx.read2contig.sam xxx.barcodeList xxx.readNameList
BIN/Sam2ReadOnContig --prefix xxx 2>>log_sam2readoncontig || exit 1
check_output xxx.read2contig

try_backup_list xxx.barcodeOnContig xxx.contigOnBarcode xxx.pe_info xxx.pe_pair
BIN/SplitInfo --prefix xxx --parse_barcode --parse_pe 2>>log_splitinfo || exit 1
check_output xxx.barcodeOnContig xxx.contigOnBarcode xxx.pe_info xxx.pe_pair
