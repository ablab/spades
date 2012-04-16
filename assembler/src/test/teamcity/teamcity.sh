#!/bin/bash

echo "### PREPROCESSING ###"

pushd ../../../
./prepare_cfg
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "prepare_cfg finished with exit code $errlvl"
    exit $errlvl
fi
popd

pushd ../../../data
./link_teamcity.sh
popd

echo "### RUNNING ###"

./$1
errlvl=$?

echo "### POSTPROCESSING ###"

pushd ../../../data
./unlink.sh
popd

echo "### TEAMCITY INVOKATION COMPLETE, EXIT CODE = $errlvl ###"

exit $errlvl
