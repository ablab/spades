#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


echo "### PREPROCESSING ###"

pushd ../../../

./prepare_cfg
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "prepare_cfg finished with exit code $errlvl"
    exit $errlvl
fi

rm bin/ -rf
./spades_compile.sh
./cpcfg
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "spades_compile finished with exit code $errlvl"
    exit $errlvl
fi

popd
pushd ../../../data
./link_teamcity.sh
popd

echo "### PRINTING CODE TO BE EXECUTED: $1 ###"

cat $1

echo "### RUNNING ###"

./$1 $2 $3
errlvl=$?

echo "### POSTPROCESSING ###"

pushd ../../../data
./unlink.sh
popd

echo "### TEAMCITY INVOKATION COMPLETE, EXIT CODE = $errlvl ###"

exit $errlvl
