#!/bin/bash

#cp src/debruijn/config.info.template src/debruijn/config.info

#./prepare_cfg
#./cpcfg
cd data/
./link_ftp.sh
cd ..
make clean
make -j 5 all
./run rd
