#!/bin/bash

set -e
cd _generic
pushd ../../../../ > /dev/null
./spades.py spades_config.info.template
popd > /dev/null

echo ""
echo "****************************************"
echo "RUNNING datasets:"
echo *
echo ""

for f in *
do
    pushd ../../../../ > /dev/null
    echo ./spades.py configs/debruijn/datasets/_generic/$f &
    ./spades.py configs/debruijn/datasets/_generic/$f > /dev/null&
    sleep 3s
    popd > /dev/null
done
