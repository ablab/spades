#!/bin/bash

string="QUAKE_CROPPED_400K";

if [ "$1" != "" ] 
then
    string="$1"
fi
echo $string
    cp ../../../data/debruijn/"$string"/K55/latest/saves/distance_estimation* .
    cp ../../../data/debruijn/"$string"/K55/latest/estimation_qual/f* .

mv fp.prd distance_estimation_fpr.prd
mv fn.prd distance_estimation_fnr.prd

