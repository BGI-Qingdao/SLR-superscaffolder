#!/bin/bash

./step_0_prepare_info.sh
./step_1_calc_seeds.sh
./step_2_bin_cluster.sh
./step_3_mst.sh
./step_4_gap_oo.sh
./step_5_trunk2scaff.sh
./clean_prepare.sh
