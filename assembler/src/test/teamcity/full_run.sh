#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
project_name=$2
creature_name=$1
project_dir="spades_output/"$creature_name$project_name
output_pref="/home/dantipov/"
options_dir=$output_pref"input/"$creature_name$project_name
rm -rf $project_dir
mkdir $project_dir -p
read line < $options_dir"spades.options"
echo $line

python ./spades.py $line -o $project_dir --disable-gzip-output 

read line < $options_dir"quast.options"
echo $line

dir=$output_pref"contigs/teamcity/"$creature_name$project_name #$(date +%Y%m%d_%H%M%S)"/"
echo $dir
mkdir $dir -p
cp $project_dir"contigs.fasta" $dir$(date +%Y%m%d_%H%M%S)".fasta"

python2.6 ~/quast/quast.py $project_dir"contigs.fasta" $line -o $project_dir"/quality_results/"
#espected results
read line < $options_dir"results.options"
echo $line
opts=( $line )
python src/test/teamcity/assess.py $project_dir"quality_results/transposed_report.tsv" ${opts[1]} ${opts[3]}

popd
