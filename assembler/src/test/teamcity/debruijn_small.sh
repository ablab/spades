#!/bin/bash
set -e
pushd ../../../
make clean
./cpcfg
./spades.py src/test/teamcity/spades_config.small.info
popd