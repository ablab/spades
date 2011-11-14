#!/bin/bash

cd ../../../
#./prepare_cfg
./cpcfg
cd data/
./link_ftp.sh
cd ..
make clean
make -j 5 all
cd src/test/teamcity/
