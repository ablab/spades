#!/bin/sh

string="ECOLI_SC_LANE_1_BH_woHUMAN";
if [ "$1" = "-p" ]
then
    path=$2
elif [ "$1" != "" ]
then
    string=$1
    path=../../../data/debruijn/$string/K55/latest
else
    path=../../../data/debruijn/$string/K55/latest
fi

#echo $path
    cp $path/saves/distance_estimation* .
    cp $path/estimation_qual/* .
    cp $path/etalon_corrected_by_graph.prd distance_estimation_et.prd 
    cp $path/etalon*.prd . 
    cp $path/scaf*.prd . 

sed '1d' distance_estimation_et.prd > etalon.prd
sed '1d' distance_estimation_cl.prd > clustered.prd

./genStats.sh

sort -rnk 4,4 fp.prd > fpr.prd
mv fpr.prd fp.prd
sort -rnk 4,4 tp.prd > tpr.prd
mv tpr.prd tp.prd
sort -rnk 4,4 etalon.prd > temp.prd
mv temp.prd etalon.prd
sort -rnk 4,4 clustered.prd > temp.prd
mv temp.prd clustered.prd

#javac PlotFPR.java
java PlotFPR -s


