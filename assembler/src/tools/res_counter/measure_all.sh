#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#### data ####

reads1=/acestorage/data/input/E.coli/is220/s_6_1.cor.fastq
reads2=/acestorage/data/input/E.coli/is220/s_6_2.cor.fastq
reads1_nocor=/acestorage/data/input/E.coli/is220/s_6_1.fastq
reads2_nocor=/acestorage/data/input/E.coli/is220/s_6_2.fastq
reads_fasta12=/acestorage/data/input/E.coli/is220/s_6.cor.fasta

## NOTE: you should add --sc option to SPAdes if your dataset is a single-cell one! ##

#### common settings ####
k=55  # where applicable
threads=16
output_base=/Johnny/tmp/measuring/E.coli/MC
res_counter=./run.pl
res_saver=./result_saver.sh


# converting to absolute paths
output_base=`readlink -e $output_base`
res_counter=`readlink -e $res_counter`
res_saver=`readlink -e $res_saver`

#### binaries #### Just comment those you don't want to start
spades_bin=/home/gurevich/misc/SPAdes-2.4.0-Linux/bin/spades.py
idba_bin=/acestorage/software/idba-1.1.0/bin/idba_ud
ray_bin=/acestorage/software/ray2/Ray
abyss_bin=abyss-pe
velvet_bin_home=/acestorage/software/velvet/
clc_bin=/acestorage/software/clc/clc_assembler
a5_bin=/acestorage/software/a5/bin/a5_pipeline.pl

#### start ####
echo "Start measuring"
mkdir -p $output_base

## SPAdes
if ! [ -z "${spades_bin+xxx}" ]; then 
echo "SPADES"
cur_out=$output_base/SPAdes2
mkdir $cur_out

$res_counter $spades_bin -1 ${reads1_nocor} -2 ${reads2_nocor} -t $threads -o $cur_out --careful
$res_saver $cur_out/measure_results
fi

## IDBA-UD
if ! [ -z "${idba_bin+xxx}" ]; then
echo "IDBA"
cur_out=$output_base/IDBA
mkdir $cur_out

$res_counter $idba_bin -r $reads_fasta12 --num_threads $threads -o $cur_out
$res_saver $cur_out/measure_results
fi

## Ray
if ! [ -z "${ray_bin+xxx}" ]; then
echo "Ray"
cur_out=$output_base/Ray
mkdir $cur_out

$res_counter mpirun -np $threads $ray_bin -p $reads1 $reads2 -o $cur_out/out -k $k
$res_saver $cur_out/measure_results
fi 

## ABySS
if ! [ -z "${abyss_bin+xxx}" ]; then
echo "ABySS"
cur_out=$output_base/ABySS
mkdir $cur_out
cur_dir=`pwd`

cd $cur_out
$res_counter $abyss_bin -j $threads k=$k name=out lib="'"pe1"'" pe1="'"$reads1 $reads2"'"
$res_saver $cur_out/measure_results
cd $cur_dir
fi

## Velvet
if ! [ -z "${velvet_bin_home+xxx}" ]; then
echo "Velvet"
cur_out=$output_base/Velvet
mkdir $cur_out

$res_counter $velvet_bin_home/velveth $cur_out $k -shortPaired -separate -fastq $reads1 $reads2
$res_counter $velvet_bin_home/velvetg $cur_out -exp_cov auto -cov_cutoff auto -scaffolding no
$res_saver $cur_out/measure_results
fi

## CLC
if ! [ -z "${clc_bin+xxx}" ]; then 
echo "CLC"
cur_out=$output_base/CLC
mkdir $cur_out

$res_counter $clc_bin -o $cur_out/contigs.fasta -p fb ss 180 260 -q -i $reads1 $reads2 --cpus $threads
$res_saver $cur_out/measure_results 
fi

## A5
if ! [ -z "${a5_bin+xxx}" ]; then
echo "A5"
cur_out=$output_base/A5
mkdir $cur_out
cur_dir=`pwd` 

cd $cur_out
$res_counter $a5_bin $reads1 $reads2 out
$res_saver $cur_out/measure_results
cd $cur_dir
fi

#### finish ####
echo "Finish measuring"
