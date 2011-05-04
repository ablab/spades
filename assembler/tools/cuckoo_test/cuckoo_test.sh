#!/bin/bash

DATA_DIR=~/Assembler/data
EXEC_DIR=../../build/tools/cuckoo_test
FILE=s_6.first10000_1.fastq.gz 

make
cp diagrams.gnu $EXEC_DIR
cd $EXEC_DIR
rm -rf temp.tmp time_insert.tmp memory.tmp
echo "File:" $DATA_DIR/$FILE 
./cuckoo_test $DATA_DIR/$FILE > temp.tmp
cat temp.tmp | grep Insert > time_insert.tmp 
cat temp.tmp | grep Memory > memory.tmp 

gnuplot 'diagrams.gnu'

echo "Diagrams are done!"
