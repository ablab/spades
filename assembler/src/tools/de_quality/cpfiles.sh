#!/bin/sh

string="ECOLI_SC_LANE_1_BH_woHUMAN";

if [ "$1" != "" ] 
then
    string="$1"
fi
echo $string
    cp ../../../data/debruijn/"$string"/K55/latest/saves/distance_estimation* .
    #cp ../../../data/debruijn/"$string"/K55/latest/estimation_qual/* .
    cp ../../../data/debruijn/"$string"/K55/latest/etalon_paired_corrected.prd distance_estimation_et.prd 

sed '1d' distance_estimation_et.prd > etalon.prd
sed '1d' distance_estimation_cl.prd > clustered.prd

./genStats.sh

sort -rnk 4,4 fp.prd > fpr.prd
mv fpr.prd fp.prd
sort -rnk 4,4 tp.prd > tpr.prd
mv tpr.prd tp.prd
sort -rnk 3,4 etalon.prd > temp.prd
mv temp.prd etalon.prd
sort -rnk 3,4 clustered.prd > temp.prd
mv temp.prd clustered.prd
