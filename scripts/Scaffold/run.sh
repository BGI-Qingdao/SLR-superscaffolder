#!/bin/bash

test -e "./__end_step1.txt" || ( echo `date`>>"./__start_step1.txt"\
  && ./step_1_prepare_info.sh && echo `date` >>"./__end_step1.txt" ) || exit 1

test -e "./__end_step2.txt" || ( echo `date`>>"./__start_step2.txt"\
   && ./step_2_order.sh && echo `date` >>"./__end_step2.txt" ) || exit 1

test -e "./__end_step3.txt" || ( echo `date`>>"./__start_step3.txt"\
   && ./step_3_gap_oo.sh && echo `date` >>"./__end_step3.txt" ) || exit 1

test -e "./__end_step4.txt" || ( echo `date`>>"./__start_step4.txt"\
   && ./step_4_pe_fill.sh && echo `date` >>"./__end_step4.txt" ) || exit 1

test -e "./__end_step5.txt" || ( echo `date`>>"./__start_step5.txt"\
   && ./step_5_gapsize.sh && echo `date` >>"./__end_step5.txt" ) || exit 1

test -e "./__end_step6.txt" || ( echo `date`>>"./__start_step6.txt"\
   && ./step_6_gen_seq.sh && echo `date` >>"./__end_step6.txt" ) || exit 1

