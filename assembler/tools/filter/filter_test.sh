#!/bin/bash

DATA_DIR=~/Assembler/data
EXEC_DIR=../../build/tools/filter
FILE_NUM=3
FILE1=s_6.first10000_1.fastq.gz 
FILE2=s_6.first100000_1.fastq.gz 
FILE3=s_6.first400000_1.fastq.gz 

make
cp diagrams.gnu $EXEC_DIR
cd $EXEC_DIR
rm -rf temp.tmp time.tmp memory.tmp
 
for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
  do echo "File:" $f; 
  echo " " > temp.tmp
  for n in 1 2 3 4 5; 
    do echo "Map type:" $n;
    (time ./filter $f 1 $n --stat >> temp.tmp) 2>>temp.tmp
  done
  cat temp.tmp | grep real | sed 's/[ms]/\t/g' >> time.tmp 
  echo " " >> time.tmp
  echo " " >> time.tmp
  cat temp.tmp | grep Memory >> memory.tmp 
  echo " " >> memory.tmp
  echo " " >> memory.tmp
done

gnuplot 'diagrams.gnu'

echo "Diagrams are done!"
