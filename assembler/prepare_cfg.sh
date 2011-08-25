#!/bin/bash

# compile and install libraries from ext

mkdir -p build/ext/lib

cd ./ext/src/io_lib-1.12.5
if [ "$1" == "nosys" ]; then
./configure --prefix="`pwd`/../../../build/ext"
else
./configure 
fi
make
make install
make clean

cd ../statgen/lib
make
cp ./libStatGen.a ../../../../build/ext/lib
cp ./samtools/libbam.a ../../../../build/ext/lib
cd ../../../include/statgen/lib/
rm -f include
ln -s ../../../src/statgen/lib/include include
cd ../samtools
rm -f bam.h
ln ../../../src/statgen/lib/samtools/bam.h bam.h
cd ../../../.. #back to assembler

# prepare debug

mkdir -p build/debug
cd build/debug
cmake ../../src -Dcfg="Debug"
cd ../.. # back to assembler

# prepare release

mkdir -p build/release
cd build/release
cmake ../../src
