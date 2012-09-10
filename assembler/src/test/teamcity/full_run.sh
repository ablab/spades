






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
./cpcfg
python ./spades.py $line -o $project_dir --disable-gzip-output 

read line < $options_dir"quast.options"
echo $line

dir=$output_pref"contigs/teamcity/"$creature_name$project_name #$(date +%Y%m%d_%H%M%S)"/"
echo $dir
mkdir $dir -p
cp $project_dir"contigs.fasta" $dir$(date +%Y%m%d_%H%M%S)".fasta"

python2.6 ~/quast-1.1/quast.py $project_dir"contigs.fasta" $line -o $project_dir"/quality_results/"
rm $dir"/quast_all" -rf
#echo $dir"*.fasta" $line -o $dir"/quast_all"
dirtmp=$dir"tmp/"
rm $dirtmp -rf
mkdir $dirtmp

for i in $dir*.fasta ; do 
#find DIR -type f | while read FILENAME; do
  flag="False"
  for j in $dirtmp*.fasta ; do
    echo $i
    echo $j
    if [ -f $j ] ;
    then
      compare_res=$(python2.6 ~/algorithmic-biology/assembler/src/tools/contig_analysis/compare_fasta.py $i $j)
      echo $compare_res
      if [ $compare_res = "True" ]; 
      then
        flag="True"
        break     
      fi	 
    fi
  done  
  if [ $flag = "False" ] ;
  then
    cp $i $dirtmp
  fi
done

quast_line="$output_pref/quast-1.1/quast.py $dirtmp* $line -o $dir/quast_all/"



echo "$quast_line"
python2.6 $quast_line
ssh -n antipov@194.85.238.21 mkdir -p "/var/www/teamcity_runs/$1$2" &

scp "$dir/quast_all/report.txt" "antipov@194.85.238.21:/var/www/teamcity_runs/$1$2/report.txt"

#espected results
read line < $options_dir"results.options"
echo $line
opts=( $line )
python src/test/teamcity/assess.py $project_dir"quality_results/transposed_report.tsv" ${opts[1]} ${opts[3]} ${opts[5]} ${opts[7]}


popd
