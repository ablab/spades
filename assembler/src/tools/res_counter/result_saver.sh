#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


dir_name="res_counter_`date +"%d.%m.%y_%T"`"
if [ $# -gt 0 ]; then
    dir_name=$1
fi

if [ -d $dir_name ]; then
    echo "Error: Directory already exists!"
else 
    mkdir -p $dir_name
    mv rc_stats.txt $dir_name
    #mv rc_stats.err $dir_name
    mv rc_stdout.txt $dir_name
    mv rc_stderr.txt $dir_name
    echo "Resource counter's files moved to $dir_name"
fi

