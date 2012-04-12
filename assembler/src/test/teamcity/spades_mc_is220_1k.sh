#!/bin/bash
set -e
./prepare.sh
pushd ../../../
./spades.py src/test/teamcity/spades_config_mc_is220_1k.info
popd