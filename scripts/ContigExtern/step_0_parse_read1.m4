#!/bin/bash

STEP="step_0 "

source ./__common_function.sh || exit 1


check_input R1 R2
try_backup_list xxx.barcodeList  xxx.readNameList
BIN/ParseReadName  --read1  R1 --prefix  xxx 2>>log_parse_read_name ||exit 1
check_output  xxx.barcodeList  xxx.readNameList
