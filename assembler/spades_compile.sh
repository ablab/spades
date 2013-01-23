#!/bin/bash

PREFIX=`pwd`
BUILD_DIR=build_spades

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
set -e
cd $BUILD_DIR
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=$PREFIX ../src $*
make -j 8
make install
cd $PREFIX
