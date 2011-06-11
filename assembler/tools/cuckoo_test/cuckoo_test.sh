#!/bin/bash

DATA_DIR=../../../data/input
EXEC_DIR=../../build/tools/cuckoo_test

#FILE1=s_6.first1000_1.fastq.gz 
#FILE2=s_6.first1000_1.fastq.gz 
#FILE3=s_6.first1000_1.fastq.gz 
FILE1=s_6.first10000_1.fastq.gz 
FILE2=s_6.first100000_1.fastq.gz 
FILE3=s_6.first400000_1.fastq.gz 

COM="MAP1='2';MAP2='3';MAP3='4';MAP4='5';MAP5='6';"

make
cp diagrams.gnu $EXEC_DIR
cd $EXEC_DIR
rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp

step=12
for d in 2 3 4 5 6; 
  do echo "Number of arrays (d):" $d;
  echo " " > temp.tmp
  for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
    do echo "File:" $f; 
    ./cuckoo_test $f $d $step ${f##*/} >> temp.tmp
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

#echo "File:" $DATA_DIR/$FILE 
#./cuckoo_test $DATA_DIR/$FILE > temp.tmp
#cat temp.tmp | grep Insert > time_insert.tmp 
#cat temp.tmp | grep Memory > memory.tmp 
#gnuplot 'diagrams.gnu'

echo "Diagrams are done!"
