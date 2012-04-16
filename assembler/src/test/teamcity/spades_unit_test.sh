#!/bin/bash
set -e
pushd ../../../
./gen_k 55
make dt
./run dt
popd
