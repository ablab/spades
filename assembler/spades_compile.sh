#!/bin/sh

############################################################################
# Copyright (c) 2015-2020 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

set -e

if [ "x$PREFIX" = "x" ]; then
  PREFIX="$(pwd)"
fi

# [misprint detecting] uncomment to check whether all of the used variables is initialized
# set -u

# default values
DEBUG="n"
RUN_TESTS="n"
BUILD_INTERNAL="n"
AMOUNT_OF_THREADS="8"
REMOVE_BUILD_DIR_BEFORE_BUILDING="n"
ADDITIONAL_FLAGS=
BUILD_CONFIG_NAME="build_config"

## there are some helpful functions ##

# https://stackoverflow.com/a/16444570
check_whether_OPTARG_is_an_integer() {
  case $OPTARG in
    (*[!0-9]*|'') echo "The argument for -j option must be an integer"; exit 2;;
    (*)           ;; # an integer
  esac
}

# return the argument first character
str_head() { echo "$(expr substr "$1" 1 1)"; }

# return the argument without the first character
str_tail() { echo "$(expr substr "$1" 2 $((${#1}-1)))"; }

print_help() {
  echo
  echo "usage:"
  echo "  $0 [options]"
  echo
  echo "options:"
  echo "  -a, --all       build internal projects"
  echo "  -d, --debug     use debug build"
  echo "  -t, --test      run basic tests (with -a runs all tests)"
  echo "  -r, --remove    remove the build directory before rebuilding"
  echo "  -j <int>        amount of threads"
  echo "  -h, --help      print this help"
  echo
  echo "  Any other options will be passed to the cmake!"
  echo
  echo "examples:"
  echo "  - build SPAdes with internal projects and run all tests, use 15 threads"
  echo "    $0 -atj15"
  echo
  echo "  - build SPAdes and run basic tests, use 9 threads"
  echo "    $0 -j 9 -t"
}

opt_a() { BUILD_INTERNAL="y"; }
opt_d() { DEBUG="y"; }
opt_t() { RUN_TESTS="y"; }
opt_r() { REMOVE_BUILD_DIR_BEFORE_BUILDING="y"; }
opt_h() { print_help; exit 0; }
opt_j() { check_whether_OPTARG_is_an_integer; AMOUNT_OF_THREADS=$OPTARG; skip_next_iter="y"; }
unknown_opt() { ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS $opt"; }

parse_long_arg() {
  case "$1" in
    (all)    opt_a;;
    (debug)  opt_d;;
    (test)   opt_t;;
    (remove) opt_r;;
    (help)   opt_h;;
    (*)      unknown_opt;; # unknown long argument
  esac
}

# returns "ok" iff the arg contains only allowed characters
verify_short_args() {
  if [ -z "$@" ]; then
    # here the $1 length is zero
    echo "ok"
  else
    # here the $1 length is non-zero
    case "$(str_head "$1")" in
      (a) echo "$(verify_short_args "$(str_tail "$1")")";;
      (d) echo "$(verify_short_args "$(str_tail "$1")")";;
      (t) echo "$(verify_short_args "$(str_tail "$1")")";;
      (r) echo "$(verify_short_args "$(str_tail "$1")")";;
      (h) echo "$(verify_short_args "$(str_tail "$1")")";;
      (j) echo "ok";; # after 'j' there should be only <int>
      (*) echo "fail";;
    esac
  fi
}

parse_short_args_handler() {
  if [ -n "$@" ]; then
    # here the $1 length is non-zero
    case "$(str_head "$1")" in
      (a) opt_a; parse_short_args_handler "$(str_tail "$1")";;
      (d) opt_d; parse_short_args_handler "$(str_tail "$1")";;
      (t) opt_t; parse_short_args_handler "$(str_tail "$1")";;
      (r) opt_r; parse_short_args_handler "$(str_tail "$1")";;
      (h) opt_h; parse_short_args_handler "$(str_tail "$1")";;
      (j) OPTARG="$(str_tail "$1")"
          if [ -z "$OPTARG" ]; then OPTARG="$next_opt"; fi
          opt_j;;
      (*) echo "WTF???"; exit 3;;
    esac
  fi
}

parse_short_args() {
  if [ "$(verify_short_args "$1")" = "ok" ]; then
    parse_short_args_handler "$1"
  else
    unknown_opt
  fi
}

get_current_build_params() {
  echo "$DEBUG $BASEDIR $ADDITIONAL_FLAGS"
}

read_buid_params() {
  if [ -z "$1" ]; then
    # there is no config folder path
    exit 5
  fi
  config_path="$1"/"$BUILD_CONFIG_NAME"
  if [ -e "$config_path" ]; then
    echo "$(cat "$config_path")"
  else
    echo ""
  fi
}

save_build_params() {
  if [ -z "$1" ]; then
    # there is no config folder path
    exit 5
  fi
  config_path="$1"/"$BUILD_CONFIG_NAME"
  echo "$(get_current_build_params)" > "$config_path"
}

normalize_path() {
  current_path="$(pwd)"
  cd "$1"
  normalized_path="$(pwd)"
  cd "$current_path"
  echo "$normalized_path"
}

## the script's logic begins here ##

skip_next_iter="n"
opt=
for next_opt in "$@" ""; do
  if [ $skip_next_iter != "n" ]; then
    skip_next_iter="n"
    opt="$next_opt"
    continue
  fi
  if [ -n "$opt" ]; then
    if [ "$(str_head "$opt")" != "-" ]; then
      unknown_opt
    else
      opt_tail=$(str_tail "$opt")
      if [ "$(str_head "$opt_tail")" = "-" ]; then
        parse_long_arg "$(str_tail "$opt_tail")"
      else
        parse_short_args "$opt_tail"
      fi
    fi
  fi
  opt="$next_opt"
done

BUILD_DIR=build_spades
BASEDIR="$(normalize_path "$(pwd)"/"$(dirname "$0")")"

if [ $DEBUG = "y" ]; then
  BUILD_DIR="${BUILD_DIR}_debug"
fi

WORK_DIR="$BASEDIR/$BUILD_DIR"
mkdir -p "$WORK_DIR"
set -e

last_build_params="$(read_buid_params "$WORK_DIR")"
if [ "$(get_current_build_params)" != "$last_build_params" ]; then
  REMOVE_BUILD_DIR_BEFORE_BUILDING="y"
fi

if [ $REMOVE_BUILD_DIR_BEFORE_BUILDING = "y" ]; then
  # we can't remove WORK_DIR itself, because it might be a symbolic link
  # and we should not remove any ".*" files, because we didn't create them
  rm -rf "${WORK_DIR:?}/"*
fi

save_build_params "$WORK_DIR"

if [ $DEBUG = "y" ]; then
  ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS -DCMAKE_BUILD_TYPE=Debug"
fi

if [ $BUILD_INTERNAL = "y" ]; then
  ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS -DSPADES_BUILD_INTERNAL=ON"
else
  ADDITIONAL_FLAGS="$ADDITIONAL_FLAGS -DSPADES_BUILD_INTERNAL=OFF"
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
