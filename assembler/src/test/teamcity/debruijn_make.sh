#!/bin/bash

pushd ../../../
make clean
#./prepare_cfg
./cpcfg
pushd data/
./link_morality.sh
popd
make -j 5 rd
popd
