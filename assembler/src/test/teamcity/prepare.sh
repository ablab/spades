#!/bin/bash
set -e
pushd ../../../
./prepare_cfg
pushd data
./link_teamcity.sh
popd
popd
