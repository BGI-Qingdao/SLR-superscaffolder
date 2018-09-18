#!/bin/bash

bwa index xxx.contig >log_bwa_index 2>&1

bwa mem -t THREADS -k eval(TAIL+1) xxx.contig R1 R2 >xxx.read2contig.sam 2>log_bwa_mem

BIN/Sam2ReadOnContig --prefix xxx 2>log_sam2readoncontig

BIN/SplitInfo --prefix xxx --parse_barcode 2>log_splitinfo
