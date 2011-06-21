#!/bin/bash

DATA_DIR=../../../data/input
EXEC_DIR=../../build/tools/cuckoo_test

#FILE1=s_6.first1000_1.fastq.gz 
#FILE2=s_6.first1000_1.fastq.gz 
#FILE3=s_6.first1000_1.fastq.gz 
FILE1=s_6.first10000_1.fastq.gz 
FILE2=s_6.first100000_1.fastq.gz 
FILE3=s_6.first400000_1.fastq.gz 

make
cp diagrams.gnu $EXEC_DIR
cd $EXEC_DIR

rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp
lbs=0
step=12
mlf=2
for d in 2 3 4 5 6; 
  do echo "Number of arrays (d):" $d;
  echo " " > temp.tmp
  for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
    do echo "File:" $f; 
    ./cuckoo_test $f $d $lbs $step $mlf ${f##*/} >> temp.tmp
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

COM="MAP1='2';MAP2='3';MAP3='4';MAP4='5';MAP5='6';FN_INS='d_t_ins.png';FN_FIND='d_t_find.png';FN_MEM='d_mem.png';"
gnuplot -e $COM diagrams.gnu

rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp
d=3
lbs=0
mlf=2
for step in 11 12 15 18 20; 
  do echo "step*10:" $step;
  echo " " > temp.tmp
  for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
    do echo "File:" $f; 
    ./cuckoo_test $f $d $lbs $step $mlf ${f##*/} >> temp.tmp
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

COM="MAP1='1.1';MAP2='1.2';MAP3='1.5';MAP4='1.8';MAP5='2.0';FN_INS='step_t_ins.png';FN_FIND='step_t_find.png';FN_MEM='step_mem.png';"
gnuplot -e $COM diagrams.gnu

rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp
d=3
lbs=0
step=15
for mlf in 1 2 5 10 20; 
  do echo "max_loop_factor:" $mlf;
  echo " " > temp.tmp
  for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
    do echo "File:" $f; 
    ./cuckoo_test $f $d $lbs $step $mlf ${f##*/} >> temp.tmp
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

COM="MAP1='1';MAP2='2';MAP3='5';MAP4='10';MAP5='20';FN_INS='mlf_t_ins.png';FN_FIND='mlf_t_find.png';FN_MEM='mlf_mem.png';"
gnuplot -e $COM diagrams.gnu

step=15
mlf=10
for d in 2 3 4;
  do rm -rf temp.tmp time_insert.tmp time_find.tmp memory.tmp
  echo "d:" $d;
  for lbs in 0 1 2 3; 
    do echo "log_2(bucket_size):" $lbs;
    echo " " > temp.tmp
    for f in $DATA_DIR/$FILE1 $DATA_DIR/$FILE2 $DATA_DIR/$FILE3; 
      do echo "File:" $f; 
      ./cuckoo_test $f $d $lbs $step $mlf ${f##*/} >> temp.tmp
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
COM="MAP1='0';MAP2='1';MAP3='2';MAP4='3';MAP5='4';FN_INS='lbs_"$d"_t_ins.png';FN_FIND='lbs_"$d"_t_find.png';FN_MEM='lbs_"$d$"_mem.png';"
gnuplot -e $COM diagrams.gnu
done

echo "Diagrams are done!"
