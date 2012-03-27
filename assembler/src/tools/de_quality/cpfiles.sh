#!/bin/sh

string="SC_HAMMER_1410_FILTERED";

if [ "$1" != "" ] 
then
    string="$1"
fi
echo $string
    cp ../../../data/debruijn/"$string"/K55/latest/saves/distance_estimation* .
    cp ../../../data/debruijn/"$string"/K55/latest/estimation_qual/* .
    cp ../../../data/debruijn/"$string"/K55/latest/etalon_paired_corrected.prd distance_estimation_et.prd 

cp fp.prd distance_estimation_fpr.prd
cp fn.prd distance_estimation_fnr.prd

sed '1d' fp.prd > temp.prd
mv temp.prd fp.prd
sed '1d' pm.prd > temp.prd
mv temp.prd pm.prd
sort -rnk 4,4 fp.prd > fpr.prd
sort -rnk 4,4 pm.prd > tp.prd
sort -rnk 3,4 distance_estimation_et.prd > etalon.prd
sort -rnk 3,4 distance_estimation_cl.prd > clustered.prd
