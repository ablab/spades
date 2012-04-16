#!/bin/bash

set -e
./../../../src/tools/spades_pipeline/generate_config.py _generic.template $*
cd _generic
for f in *
do
    pushd ../../../../
    echo ./spades.py src/debruijn/datasets/_generic/$f &
    #./spades.py src/debruijn/datasets/_generic/$f &
    popd
done