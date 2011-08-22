#!/bin/bash

mkdir -p ./build/libs

cd ./ext/src/io_lib-1.12.5
./configure --prefix=../../../build/libs
make
make install
make clean

cd ../statgen/lib
make
cp ./libStatGen.a ../../../../build/libs
cp ./samtools/libbam.a ../../../../build/libs
