#!/bin/bash
#

test -e "./__end_step0.txt" || ( echo `date`>>"./__start_step0.txt"\
  && ./step_0_parse_read1.sh && echo `date` >>"./__end_step0.txt" ) || exit 1

test -e "./__end_step1.txt" || ( echo `date`>>"./__start_step1.txt"\
   && ./step_1_prepare_info.sh && echo `date` >>"./__end_step1.txt" ) || exit 1

test -e "./__end_step2.txt" || ( echo `date`>>"./__start_step2.txt"\
   && ./step_2_calc_seeds.sh && echo `date` >>"./__end_step2.txt" ) || exit 1

test -e "./__end_step3.txt" || ( echo `date`>>"./__start_step3.txt"\
   && ./step_3_extern_contig.sh && echo `date` >>"./__end_step3.txt" ) || exit 1

test -e "./__end_step4.txt" || ( echo `date`>>"./__start_step4.txt"\
   && ./step_4_merge_contig.sh && echo `date` >>"./__end_step4.txt" ) || exit 1

