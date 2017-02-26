#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e

if [ "x$PREFIX" = "x" ]; then
  PREFIX=`pwd`
fi
BUILD_DIR=build_spades
BASEDIR=`pwd`/`dirname $0`

rm -rf "$BASEDIR/$BUILD_DIR"
mkdir -p "$BASEDIR/$BUILD_DIR"
set -e
cd "$BASEDIR/$BUILD_DIR"
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$PREFIX" $* "$BASEDIR/src"
make -j 8
make install
cd $PREFIX
