#!/bin/bash 
set -e
./prepare.sh
pushd ../../../
./gen_k 55
make dt
./run dt
popd
