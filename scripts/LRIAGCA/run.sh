#!/bin/bash
#

test -e "./__end_step0.txt" || ( echo `date`>>"./__start_step0.txt"\
  && ./step_0_prepare_info.sh && echo `date` >>"./__end_step0.txt" ) || exit 1

test -e "./__end_step1.txt" || ( echo `date`>>"./__start_step1.txt"\
   && ./step_1_calc_seeds.sh && echo `date` >>"./__end_step1.txt" ) || exit 1

test -e "./__end_step2.txt" || ( echo `date`>>"./__start_step2.txt"\
   && ./step_2_bin_cluster.sh && echo `date` >>"./__end_step2.txt" ) || exit 1

test -e "./__end_step3.txt" || ( echo `date`>>"./__start_step3.txt"\
   && ./step_3_contig_dlink.sh && echo `date` >>"./__end_step3.txt" ) || exit 1

test -e "./__end_step4.txt" || ( echo `date`>>"./__start_step4.txt"\
   && ./step_4_contig_road.sh && echo `date` >>"./__end_step4.txt" ) || exit 1

test -e "./__end_step5.txt" || ( echo `date`>>"./__start_step5.txt"\
   && ./step_5_contig_road_fill.sh && echo `date` >>"./__end_step5.txt" ) || exit 1

test -e "./__end_step6.txt" || ( echo `date`>>"./__start_step6.txt"\
   && ./step_6_merge_contig.sh && echo `date` >>"./__end_step6.txt" ) || exit 1
