#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "script.sh <in_prefix> <out_dir>"
    exit
fi

in_prefix=$1
out_dir=$2

for d in $in_prefix* ; do
   name=$(basename $d) 
   fn=$out_dir/$name.desc
   echo $d/left.fastq > $fn
   echo $d/right.fastq >> $fn
done
