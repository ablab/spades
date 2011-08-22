#!/bin/bash

mkdir -p ./build/libs

cd ./ext/src/io_lib-1.12.5
./configure
make
make install
make clean

cd ../statgen/lib
make
cp ./libStatGen.a ../../../../build/libs
cp ./samtools/libbam.a ../../../../build/libs
cd ../../../include/statgen/lib/
rm -f include
ln -s ../../../src/statgen/lib/include include
cd ../samtools
rm -f bam.h
ln ../../../src/statgen/lib/samtools/bam.h bam.h