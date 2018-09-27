STEP="step_5 "

source ./__common_function.sh || exit 1

check_input xxx.contig
try_backup_list xxx.seeds
BIN/StaticsticUnique --prefix xxx --kvalue KVALUE --min PE_SEED_MIN 2>>log_staticsticunique
check_output xxx.seeds

check_input xxx.pe_pairs xxx.pe_info 
try_backup_list xxx.pe_graph
BIN/PEGraph --prefix xxx --max_is MAX_IS  2>>log_pe_graph
check_output xxx.pe_graph

check_input xxx.seeds xxx.barcodeOnContig
try_backup_list xxx.barcodeOnBin
BIN/ChopBin --prefix xxx --bin_size PE_BIN_SIZE  --delete_tail TAIL  2>>log_chopbin_pe
check_output xxx.barcodeOnBin

check_input xxx.barcodeOnBin
try_backup_list xxx.cluster
BIN/BinCluster --prefix xxx --thread THREADS --threshold PE_BIN_THRESHOLD 2>>log_bincluster_pe
check_output xxx.cluster

check_input xxx.mintree_trunk_linear xxx.cluster
try_backup_list xxx.seeds_cluster_seeds
BIN/SeedCluster --prefix xxx  --threshold PE_BIN_THRESHOLD --strategy 1  --loop_num 1 2>>log_seed_cluster
check_output xxx.seeds_cluster_seeds

check_input xxx.seeds xxx.seeds_cluster_seeds xxx.pe_graph
try_backup_list xxx.trunk_fill
BIN/FillTrunkByPE  --prefix xxx  --searchMax  PE_SEARCH_MAX --insert_max MAX_IS  --min_count  PE_MIN_COUNT 2>>log_fill_trunk_by_pe 
check_output xxx.trunk_fill
