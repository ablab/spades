#!/bin/sh

SEQTKPATH=/home/tdvorkina/soft/seqtk/

python2 filterbylength.py $1 ./benchmark_data/ecoli/tmp/${2}2000_all.$3 2000
$SEQTKPATH/seqtk sample -s100 ./benchmark_data/ecoli/tmp/${2}2000_all.$3 10000 > ./benchmark_data/ecoli/input/${2}2000.$3
rm -rf ./benchmark_data/ecoli/tmp/
