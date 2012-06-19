#!/bin/bash

## WARNING! Should be rewritten ASAP! Very very strange paired mode

BUILD_DIR=${PWD}/../../../build/ext/velvet

if [ ! -d ${BUILD_DIR} ]
then 
    echo "== copying Velvet to build dir =="
    mkdir -p ${BUILD_DIR}
    cp -r * ${BUILD_DIR}
    echo "== copying finished =="
fi

cd ${BUILD_DIR}

if [ ! -f velvetg ] || [ ! -f velveth ]
then 
    echo "== making Velvet =="
    make >${BUILD_DIR}/make.log 2>${BUILD_DIR}/make.err
    echo "== making finished =="
fi

echo "== assembling started =="
OUTPUT_DIR=${BUILD_DIR}/results
rm -rf $OUTPUT_DIR
mkdir tmp

for i in $*
do
  cp $i tmp
done
gzip -d tmp/*
mkdir tmp/mixed
if [ -e "tmp/*left.corrected*.fastq" ]
then
  perl shuffleSequences_fastq.pl tmp/*left.corrected*.fastq tmp/*right.corrected*.fastq tmp/mixed/mixed.fastq
  rm tmp/*corrected*.fastq
elif [ -e "tmp/*left.corrected*.fasta" ]
then
  perl shuffleSequences_fasta.pl tmp/*left.corrected*.fasta tmp/*right.corrected*.fastq tmp/mixed/mixed.fasta
  rm tmp/*corrected*.fasta

elif [ -e "tmp/*_1.fastq" ]
then
  echo "here"
  perl shuffleSequences_fastq.pl tmp/*_1.fastq tmp/*_2.fastq tmp/mixed/mixed.fastq
  rm tmp/*_1.fastq
  rm tmp/*_2.fastq
elif [ -e "tmp/*_1.fasta" ]
then
  perl shuffleSequences_fasta.pl tmp/*_1.fasta tmp/*_2.fasta tmp/mixed/mixed.fasta
  rm tmp/*_1.fasta
  rm tmp/*_2.fasta
fi
velvethrun="velveth results 55"
if [ -e "tmp/mixed/mixed.fasta" ]
then
  let "velvethrun += " -shortPaired -fasta tmp/mixed/mixed.fasta""
elif [ -e "tmp/mixed/mixed.fastq" ]
then
  let "velvethrun += " -shortPaired -fastq tmp/mixed/mixed.fastq""
fi
#for filename in tmp/;
FOUND=0

echo $velvethrun
cd tmp/

for filename in *.f*a        # Перебор всех файлов в текущем каталоге. (fasta, fa)
do
echo $filename
if [ -e $filename ]
then
velvethrun="$velvethrun -short -fasta tmp/$filename"
fi
done
for filename in *.f*q        # Перебор всех файлов в текущем каталоге. (fastq, fq)
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
./velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -scaffolding no
rm tmp -rf

#perl shuffleSequences_fastq.pl tmp/*_1*.fastq tmp/*_2*.fastq tmp/mixed.fastq

echo "== assembling finished. Contigs are $OUTPUT_DIR/contigs.fa =="
