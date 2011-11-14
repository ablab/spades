#!/bin/bash

./debruijn_make.sh
cd ../../../
cp ./src/test/teamcity/config.info.FULL ./src/debruijn/config.info
./run rd
