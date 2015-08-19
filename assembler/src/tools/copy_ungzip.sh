#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

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