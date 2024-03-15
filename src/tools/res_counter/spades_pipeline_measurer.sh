#!/bin/bash

############################################################################
# Copyright (c) 2015-2016 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# reads options (you probably need to change it!)
sc=""  # left empty if multicell
reads_base_dir=/Nancy/data/input/Bacteria/E.coli/K12/is220
left=$reads_base_dir/s_6_1.fastq.gz
right=$reads_base_dir/s_6_2.fastq.gz

# reads options (you probably need to change it!)
#sc="--sc"  # left empty if multicell
#reads_base_dir=/Nancy/data/input/Bacteria/E.coli/K12/ucsd_lane_1
#left=$reads_base_dir/ecoli_mda_lane1_left.fastq
#right=$reads_base_dir/ecoli_mda_lane1_right.fastq


# or you can set reads via command line arguments ("left right --sc" or "left right")
if [[ -f $1 && -f $2 ]]; then
  left=$1
  right=$2
  sc=$3
fi

# advanced options (if you are at @rep/algorithmic-biology/assembler/, these params should work for you)
spades_dir=.  # use . if your cwd contains spades.py
res_counter_dir=./src/tools/res_counter
save_results_here=./stats_results
spades_output_dir=./spades_output
other_spades_py_params="-k 21,33,55 --careful"

# for memchecker, don't change if you are not sure!
stopper_fname=mem_checker_stopper.txt  # should be the same to fname in mem_checker.py
timeout=7  # in seconds, should be >= timeout in mem_checker.py

echo "STARTING, final stats will be here: $save_results_here"
rm -fr $spades_output_dir
touch $stopper_fname
echo "ERROR CORRECTION"
$res_counter_dir/mem_checker.py $spades_output_dir &
$res_counter_dir/run.pl $spades_dir/spades.py -o $spades_output_dir $other_spades_py_params $sc -1 $left -2 $right --stop-after ec
rm $stopper_fname
sleep $timeout
$res_counter_dir/result_saver.sh $save_results_here/ec_stats
echo -n "disk peak usage: " >> $save_results_here/ec_stats/rc_stats.txt
cat maxmem.txt >> $save_results_here/ec_stats/rc_stats.txt
rm maxmem.txt
previous_stage_size=`du -hs $spades_output_dir | cut -f1`

touch $stopper_fname
echo "ASSEMBLING"
$res_counter_dir/mem_checker.py $spades_output_dir  &
$res_counter_dir/run.pl $spades_dir/spades.py -o $spades_output_dir --restart-from last --stop-after as
rm $stopper_fname
sleep $timeout
$res_counter_dir/result_saver.sh $save_results_here/as_stats
echo -n "disk peak usage: " >> $save_results_here/as_stats/rc_stats.txt
echo "Don't forget to exclude disk usage by previous stage results: $previous_stage_size!" >> maxmem.txt
cat maxmem.txt >> $save_results_here/as_stats/rc_stats.txt
rm maxmem.txt
previous_stage_size=`du -hs $spades_output_dir | cut -f1`

touch $stopper_fname
echo "MISMATCH CORRECTOR"
$res_counter_dir/mem_checker.py $spades_output_dir  &
$res_counter_dir/run.pl $spades_dir/spades.py -o $spades_output_dir --restart-from last --stop-after mc
rm $stopper_fname
sleep $timeout
$res_counter_dir/result_saver.sh $save_results_here/mc_stats
echo -n "disk peak usage: " >> $save_results_here/mc_stats/rc_stats.txt
echo "Don't forget to exclude disk usage by previous stage results: $previous_stage_size!" >> maxmem.txt
cat maxmem.txt >> $save_results_here/mc_stats/rc_stats.txt
rm maxmem.txt

final_stats=$save_results_here/stats.txt
echo "Pipeline stats: " > $final_stats
echo "" >> $final_stats
echo "Error correction: " >> $final_stats
cat $save_results_here/ec_stats/rc_stats.txt >> $final_stats
echo "" >> $final_stats
echo "Assembling: " >> $final_stats
cat $save_results_here/as_stats/rc_stats.txt >> $final_stats
echo "" >> $final_stats
echo "Mismatch corrector: " >> $final_stats
cat $save_results_here/mc_stats/rc_stats.txt >> $final_stats
echo "DONE! All stats are concatenated here: $final_stats"
