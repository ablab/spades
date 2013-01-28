#!/bin/bash

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

PREFIX=`pwd`
BUILD_DIR=build_spades
BIN_DIR=bin

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
set -e
cd $BUILD_DIR
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=$PREFIX ../src $*
make -j 8
make install
cd $PREFIX

# bwa for mismatch corrector
BWA_SRC_DIR=./ext/tools/bwa-0.6.2
BWA_BUILD_DIR=$BUILD_DIR/bwa
mkdir -p $BUILD_DIR
cp -r $BWA_SRC_DIR $BWA_BUILD_DIR
echo -e "\nBuilding BWA for MismatchCorrector tool"
make -C $BWA_BUILD_DIR
cp $BWA_BUILD_DIR/bwa $BIN_DIR
