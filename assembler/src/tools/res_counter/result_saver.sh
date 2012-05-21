#!/bin/bash

dir_name="res_counter_`date +"%d.%m.%y_%T"`"
if [ $# -gt 0 ]; then
    dir_name=$1
fi

if [ -d $dir_name ]; then
    echo "Error: Directory already exists!"
else 
    mkdir -p $dir_name
    mv rc_stats.txt $dir_name
    mv rc_stats.err $dir_name
    mv rc_stdout.txt $dir_name
    mv rc_stderr.txt $dir_name
	echo "Resource counter's files moved to $dir_name"
fi

