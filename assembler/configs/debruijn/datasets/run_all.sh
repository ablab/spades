#!/bin/bash

set -e
cd _generic
for f in *
do
    pushd ../../../../ > /dev/null
    echo ./spades.py configs/debruijn/datasets/_generic/$f &
    ./spades.py configs/debruijn/datasets/_generic/$f > /dev/null&
    sleep 5s
    popd > /dev/null
done
