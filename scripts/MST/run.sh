#!/bin/bash

test -e "./__end_step1.txt" || ( echo `date`>>"./__start_step1.txt"\
  && ./step_1_prepare_info.sh && echo `date` >>"./__end_step1.txt" ) || exit 1

test -e "./__end_step2.txt" || ( echo `date`>>"./__start_step2.txt"\
   && ./step_2_bin_cluster.sh && echo `date` >>"./__end_step2.txt" ) || exit 1

test -e "./__end_step3.txt" || ( echo `date`>>"./__start_step3.txt"\
   && ./step_3_oo.sh && echo `date` >>"./__end_step3.txt" ) || exit 1

test -e "./__end_step5.txt" || ( echo `date`>>"./__start_step5.txt"\
   && ./step_5_trunk2scaff.sh && echo `date` >>"./__end_step5.txt" ) || exit 1

