#!/bin/bash

# compile and install libraries from ext

export assembler="`pwd`"
export src="$assembler/src"
export ext="$assembler/ext"
export build="$assembler/build"

echo "/**************************************/"
echo "/*****  building extern libraries *****/"
echo "/**************************************/"
echo "----"

ext/prepare_ext.sh

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


echo "/**************************************/"
echo "/********  generating k.hpp  **********/"
echo "/**************************************/"
./gen_k.sh 55
