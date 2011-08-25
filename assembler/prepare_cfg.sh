#!/bin/bash

# compile and install libraries from ext

export assembler="`pwd`"
export src="$assembler/src"
export ext="$assembler/ext"
export build="$assembler/build"



#  echo "/**************************************/"
#  echo "/********  building staden ************/"
#  echo "/**************************************/"

#  mkdir -p $build/ext/staden
#  cd $build/ext/staden

#  $ext/src/io_lib-1.12.5/configure --prefix="`pwd`"

#  make
#  make install

#  echo "/**************************************/"
#  echo "/********  building statgen ***********/"
#  echo "/**************************************/"

#  mkdir -p $build/ext/statgen

#  cd $ext/src/statgen/lib
#  make

#  cp ./libStatGen.a       $build/ext/statgen/
#  cp ./samtools/libbam.a  $build/ext/statgen/

#unfortunately need to clean everything 
#  make clean 

#cd ../../../include/statgen/lib/
#rm -f include
#ln -s ../../../src/statgen/lib/include include
#cd ../samtools
#rm -f bam.h
#ln ../../../src/statgen/lib/samtools/bam.h bam.h
#cd ../../../.. #back to assembler

echo "/**************************************/"
echo "/********  preparing debug ************/"
echo "/**************************************/"

mkdir -p $build/debug
cd $build/debug
cmake $src -Dcfg="Debug"

echo "/**************************************/"
echo "/********  preparing release **********/"
echo "/**************************************/"

mkdir -p $build/release
cd $build/release
cmake $src
