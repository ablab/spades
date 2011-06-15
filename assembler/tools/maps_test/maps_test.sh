#!/bin/bash

DATA_DIR=../../../data/input
EXEC_DIR=../../build/tools/maps_test

#FILE1=s_6.first10000_1.fastq.gz 
#FILE2=s_6.first10000_1.fastq.gz 
#FILE3=s_6.first10000_1.fastq.gz 
FILE1=s_6.first10000_1.fastq.gz 
FILE2=s_6.first100000_1.fastq.gz 
FILE3=s_6.first400000_1.fastq.gz 

COM="MAP1='map';MAP2='hash_map';MAP3='unordered_map';MAP4='sparse_hash_map';MAP5='cuckoo';"

make
cp diagrams.gnu $EXEC_DIR
cd $EXEC_DIR
rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp
 
for n in 1 2 3 4 5; 
  do echo "Map type:" $n;
  echo " " > temp.tmp
  for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
    do echo "File:" $f; 
    ./filter $f $n ${f##*/} >> temp.tmp
  done
  cat temp.tmp | grep Insert >> time_insert.tmp 
  echo " " >> time_insert.tmp
  echo " " >> time_insert.tmp
  cat temp.tmp | grep Find >> time_find.tmp 
  echo " " >> time_find.tmp
  echo " " >> time_find.tmp
  cat temp.tmp | grep Memory >> memory.tmp 
  echo " " >> memory.tmp
  echo " " >> memory.tmp
done

gnuplot -e $COM diagrams.gnu

echo "Diagrams are done!"
