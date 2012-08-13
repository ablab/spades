#!/bin/bash

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e
pushd ../../../
project_dir="spades_output/ECOLI_1K/"
options_dir="/storage/data/input/E.coli.K12/1K/"
rm -rf $project_dir
read line < $options_dir"spades.options"
echo $line

python ./spades.py $line -o $project_dir
read line < $options_dir"quast.options"
echo $line

dir="/storage/data/contigs/teamcity/Ecoli1K/"$(date +%Y%m%d_%H%M%S)"/"
echo $dir
mkdir $dir
cp $project_dir"contigs.fasta" $dir

python /storage/quast1.0/quast/quast.py $project_dir"contigs.fasta" $line -o $project_dir"/quality_results/"
python src/test/teamcity/assess.py $project_dir"quality_results/transposed_report.tsv" 1000 0

popd
