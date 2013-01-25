#!/bin/bash

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

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
