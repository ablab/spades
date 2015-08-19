#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


velvethrun="velveth results 55"
if [ -e "tmp/mixed/mixed.fasta" ] 
then
  velvethrun="$velvethrun -shortPaired -fasta tmp/mixed/mixed.fasta"  
elif [ -e "tmp/mixed/mixed.fastq" ] 
then 
  velvethrun="$velvethrun -shortPaired -fastq tmp/mixed/mixed.fastq" 
fi
#for filename in tmp/; 
FOUND=0

echo $velvethrun
cd tmp/

for filename in *.fasta        # Перебор всех файлов в текущем каталоге.
do
echo $filename
if [ -e $filename ]
then
velvethrun="$velvethrun -short -fasta tmp/$filename"
fi
done
for filename in *.fastq        # Перебор всех файлов в текущем каталоге.
do
echo $filename
if [ -e $filename ]
then
velvethrun="$velvethrun -short -fastq tmp/$filename"
fi
done
cd ../
echo $velvethrun
./$velvethrun
./velvetg results -exp_cov auto -cov_cutoff auto -scaffolding no
rm tmp -rf

#perl shuffleSequences_fastq.pl tmp/*_1*.fastq tmp/*_2*.fastq tmp/mixed.fastq
