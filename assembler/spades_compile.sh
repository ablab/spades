#!/bin/bash

PREFIX=`pwd`

rm -rf build
mkdir -p build
cd build
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=$PREFIX ../src $*
make
make install
cd $PWD
