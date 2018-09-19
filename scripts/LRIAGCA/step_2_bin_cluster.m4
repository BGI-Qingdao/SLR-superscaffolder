#!/bin/bash

BIN/ChopBin --prefix xxx --bin_size BIN_SIZE  --delete_tail TAIL  2>>log_chopbin

BIN/BinCluster --prefix xxx --thread THREADS --threshold THRESHOLD 2>>log_bincluster
