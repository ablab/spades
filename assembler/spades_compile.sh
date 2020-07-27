#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2020 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e

if [ "x$PREFIX" = "x" ]; then
  PREFIX=$(pwd)
fi

# [misprint detecting] uncomment to check whether all of the used variables is initialized
# set -u

# default values
DEBUG="n"
RUN_TESTS="n"
BUILD_INTERNAL="n"
AMOUNT_OF_THREADS="8"
REMOVE_CACHE_BEFORE_BUILDING="n"
ADDITIONAL_FLAGS=

# https://stackoverflow.com/a/16444570
check_whether_OPTARG_is_an_integer() {
  case $OPTARG in
    (*[!0-9]*|'') echo "The argument for -$opt option must be an integer"; exit 2;;
    (*)           ;; # an integer
  esac
}

print_help() {
  echo
  echo "usage:"
  echo "  $0 [options]"
  echo
  echo "options:"
  echo "  -a          build internal projects"
  echo "  -d          use debug build"
  echo "  -t          run basic tests (with -a runs all tests)"
  echo "  -r          remove the build cache before rebuilding"
  echo "  -j <int>    amount of threads"
  echo
  echo "examples:"
  echo "  - build SPAdes with internal projects and run all tests, use 15 threads"
  echo "    $0 -atj15"
  echo
  echo "  - build SPAdes and run basic tests, use 9 threads"
  echo "    $0 -j 9 -t"
}

while getopts "adtrj:" opt; do
  case $opt in
    (a) BUILD_INTERNAL="y";;
    (d) DEBUG="y";;
    (t) RUN_TESTS="y";;
    (r) REMOVE_CACHE_BEFORE_BUILDING="y";;
    (j) check_whether_OPTARG_is_an_integer; AMOUNT_OF_THREADS=$OPTARG;;
    (*) print_help; exit 1;;
  esac
done

if [ $DEBUG = "y" ]; then
  ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS -DCMAKE_BUILD_TYPE=Debug"
fi

if [ $BUILD_INTERNAL = "y" ]; then
  ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS -DSPADES_BUILD_INTERNAL=ON"
fi

BUILD_DIR=build_spades
BASEDIR=$(pwd)/$(dirname "$0")

if [ $DEBUG = "y" ]; then
  BUILD_DIR="${BUILD_DIR}_debug"
fi

WORK_DIR="$BASEDIR/$BUILD_DIR"
mkdir -p "$WORK_DIR"
set -e

if [ $REMOVE_CACHE_BEFORE_BUILDING = "y" ]; then
  # we can't remove WORK_DIR itself, because it might be a symbolic link
  # and we should not remove any ".*" files, because we didn't create them
  rm -rf "${WORK_DIR:?}/"*
fi

cd "$WORK_DIR"
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$PREFIX" $ADDITIONAL_FLAGS "$BASEDIR/src"
make -j $AMOUNT_OF_THREADS
make install
cd "$PREFIX"

if [ $RUN_TESTS = "y" ]; then
  cd "$BASEDIR"
  if [ $BUILD_INTERNAL = "y" ]; then
    "$WORK_DIR/bin/include_test"   ; set -e
    "$WORK_DIR/bin/debruijn_test"  ; set -e
  fi
  SPADES="$BASEDIR"/spades.py
  "$SPADES" -t $AMOUNT_OF_THREADS --test              ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --isolate    ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --sc         ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --meta       ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --bio        ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --rna        ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --plasmid    ; set -e
  "$SPADES" -t $AMOUNT_OF_THREADS --test --iontorrent ; set -e
  cd "$PREFIX"
fi
