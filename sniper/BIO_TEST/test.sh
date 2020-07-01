#!/bin/bash
OUT=$1
CONFIG=$2
NT=$3

APP_PATH="/media/zmac/Research"

COUNT_EXE="$APP_PATH/megahit/build/megahit_core count -k 21 -m 2 --host_mem 30130241126.0 --mem_flag 1 --output_prefix $APP_PATH/megahit/build/small_out/tmp/k21/21 --num_cpu_threads ${NT} --read_lib_file $APP_PATH/megahit/build/small_out/tmp/reads.lib"
BUILD_EXE="$APP_PATH/megahit/build/megahit_core seq2sdbg --host_mem 30130241126.0 --mem_flag 1 --output_prefix $APP_PATH/megahit/build/small_out/tmp/k21/21 --num_cpu_threads ${NT} -k 21 --kmer_from 0 --input_prefix $APP_PATH/megahit/build/small_out/tmp/k21/21 --need_mercy"
ASSEMBLE_EXE="$APP_PATH/megahit/build/megahit_core assemble -s /media/zmac/SpinDisk/Research/megahit/build/small_out/tmp/k21/21 -o $APP_PATH/megahit/build/small_out/intermediate_contigs/k21 -t 8 --min_standalone 200 --prune_level 2 --merge_len 20 --merge_similar 0.95 --cleaning_rounds 5 --disconnect_ratio 0.1 --low_local_ratio 0.2 --cleaning_rounds 5
--min_depth 2 --bubble_level 2 --max_tip_len -1 --is_final_round"

echo $COUNT_EXE

../run-sniper -c $CONFIG -n $NT --roi -d $OUT -- $COUNT_EXE
#perf stat -e cache-misses $COUNT_EXE
