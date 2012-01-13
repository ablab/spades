#!/bin/bash

pushd ../../../
make clean
#./prepare_cfg
./cpcfg
./gen_k 55
pushd data/
./link_morality.sh
popd
make -j 5 rd
popd
