#!/bin/bash
make rd;
make rh;
mkdir -p bin;
rm bin/*;
cp build/release/debruijn/spades bin/spades;
cp build/release/build_hammer/hammer/hammer bin/hammer;
python spades.py "$@";

