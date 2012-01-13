#!/bin/bash 

pushd ../../../
make clean
./gen_k 55
make rdt
./run rdt
popd
