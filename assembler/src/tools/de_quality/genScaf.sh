#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


PROJ_DIR=`pwd`
OUTPUT_DIR=$PROJ_DIR/data/output/
INPUT_DIR=$PROJ_DIR/data/input/


cd $INPUT_DIR
sed '1d' scaf_etalon_*.prd > etalon.prd
sed '1d' scaf_clustered.prd > clustered.prd

cd $PROJ_DIR
./genStats.sh
cd $INPUT_DIR

#sort -rnk 4,4 fp.prd > fpr.prd
#mv fpr.prd fp.prd
#sort -rnk 4,4 tp.prd > tpr.prd
#mv tpr.prd tp.prd
#sort -rnk 4,4 etalon.prd > temp.prd
#mv temp.prd etalon.prd
#sort -rnk 4,4 clustered.prd > temp.prd
#mv temp.prd clustered.prd

cd $PROJ_DIR
./genPlot.sh


