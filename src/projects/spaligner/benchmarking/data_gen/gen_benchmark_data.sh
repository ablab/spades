#!/bin/sh

SEQTKPATH=

mkdir ../ecoli/tmp/
python2 filterbylength.py $1 ../ecoli/tmp/${2}2000_all.$3 2000
$SEQTKPATH/seqtk sample -s115249 ../ecoli/tmp/${2}2000_all.$3 10000 > ../ecoli/input/${2}2000.$3
rm -rf ../ecoli/tmp/


mkdir ../celegans/tmp/
python2 filterbylength.py $1 ../celegans/tmp/${2}2000_all.$3 2000
$SEQTKPATH/seqtk sample -s115249 ../celegans/tmp/${2}2000_all.$3 10000 > ../celegans/input/${2}2000.$3
rm -rf ../celegans/tmp/

mkdir ../scerevisiae/tmp/
python2 filterbylength.py $1 ../scerevisiae/tmp/${2}2000_all.$3 2000
$SEQTKPATH/seqtk sample -s115249 ../scerevisiae/tmp/${2}2000_all.$3 10000 > ../scerevisiae/input/${2}2000.$3
rm -rf ../scerevisiae/tmp/
