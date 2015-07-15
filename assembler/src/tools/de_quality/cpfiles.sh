#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


input_dir="data/input/"
output_dir="data/output/"
proj_dir=`pwd`
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

mkdir -p $output_dir
mkdir $input_dir > /dev/null
ant

cp $path/saves/distance_estimation* $input_dir
cp $path/estimation_qual/paths.prd $input_dir
cp $path/etalon_corrected_by_graph.prd $input_dir/distance_estimation_et.prd 
cp $path/etalon*.prd $input_dir 
cp $path/scaf*.prd $input_dir

cd $input_dir
sed '1d' distance_estimation_et.prd > etalon.prd
sed '1d' distance_estimation_cl.prd > clustered.prd
cd $proj_dir

./genFiles.sh

cd $input_dir
#sort -rnk 4,4 fp.prd > fpr.prd
#mv fpr.prd fp.prd
#sort -rnk 4,4 tp.prd > tpr.prd
#mv tpr.prd tp.prd
#sort -rnk 4,4 etalon.prd > temp.prd
#mv temp.prd etalon.prd
#sort -rnk 4,4 clustered.prd > temp.prd
#mv temp.prd clustered.prd
cd $proj_dir

./genPlot.sh
