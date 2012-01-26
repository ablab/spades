#!/bin/bash 
set -e
pushd ../../../
make clean
./gen_k 55
make dt
./run dt
popd
