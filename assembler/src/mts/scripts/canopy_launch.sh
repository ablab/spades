#!/bin/bash
/home/ygorshkov/Projects/canopy/cc.bin -n 32 -i $1 -o bin_canopy -c prof_canopy --max_canopy_dist 0.1 --max_close_dist 0.4 --max_merge_dist 0.05 --min_step_dist 0.01 --max_num_canopy_walks 3 --stop_fraction 1 --canopy_size_stats_file stat --filter_min_obs 1 --filter_max_dominant_obs 1.0

#/home/ygorshkov/Projects/canopy/cc.bin -n 32 -i canopy_mod.in -o bin_canopy -c prof_canopy --max_canopy_dist 0.1 --max_close_dist 0.4 --max_merge_dist 0.1 --min_step_dist 0.005 --max_num_canopy_walks 5 --stop_fraction 1 --canopy_size_stats_file stat
