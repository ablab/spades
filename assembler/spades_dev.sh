#!/bin/bash

set -e

make rd;
make rh;
mkdir -p bin;
rm bin/*;
cp build/release/debruijn/spades bin/spades;
cp build/release/hammer/hammer bin/hammer;
python spades.py "$@";

