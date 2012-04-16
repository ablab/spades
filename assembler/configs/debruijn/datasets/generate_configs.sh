#!/bin/bash

set -e
./../../../src/tools/spades_pipeline/generate_config.py _generic.template $*
cd _generic
for f in *
do
    pushd ../../../../ > /dev/null
    echo ./spades.py configs/debruijn/datasets/_generic/$f &
    ./spades.py configs/debruijn/datasets/_generic/$f &
    popd > /dev/null
done
