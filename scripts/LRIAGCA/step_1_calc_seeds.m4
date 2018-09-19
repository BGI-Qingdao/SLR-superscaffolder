#!/bin/bash
# 
# Choose seed contig
#

BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --min eval(BIN_SIZE+TAIL) 2>>log_staticsticunique
