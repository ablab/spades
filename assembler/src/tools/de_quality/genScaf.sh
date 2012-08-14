#!/bin/sh

sed '1d' scaf_etalon.prd > etalon.prd
sed '1d' scaf_clustered.prd > clustered.prd

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


