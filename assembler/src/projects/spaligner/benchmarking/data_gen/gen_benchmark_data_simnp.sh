#!/bin/sh

#FOR SIMULATED NANOPORE

SEQTKPATH=

mkdir ../ecoli/tmp/
python2 filter_unaligned.py $1 ../ecoli/tmp/${2}_aligned.fasta
python2 filterbylength.py ../ecoli/tmp/${2}_aligned.fasta ../ecoli/tmp/${2}2000_all.fasta 2000
$SEQTKPATH/seqtk sample -s115249 ../ecoli/tmp/${2}2000_all.fasta 10000 > ../ecoli/input/${2}2000.fasta
rm -rf ../ecoli/tmp/


mkdir ../celegans/tmp/
python2 filter_unaligned.py $1 ../celegans/tmp/${2}_aligned.fasta
python2 filterbylength.py ../celegans/tmp/${2}_aligned.fasta ../celegans/tmp/${2}2000_all.fasta 2000
$SEQTKPATH/seqtk sample -s115249 ../celegans/tmp/${2}2000_all.fasta 10000 > ../celegans/input/${2}2000.fasta
rm -rf ../celegans/tmp/


mkdir ../scerevisiae/tmp/
python2 filter_unaligned.py $1 ../scerevisiae/tmp/${2}_aligned.fasta
python2 filterbylength.py ../scerevisiae/tmp/${2}_aligned.fasta ../scerevisiae/tmp/${2}2000_all.fasta 2000
$SEQTKPATH/seqtk sample -s115249 ../scerevisiae/tmp/${2}2000_all.fasta 10000 > ../scerevisiae/input/${2}2000.fasta
rm -rf ../scerevisiae/tmp/
