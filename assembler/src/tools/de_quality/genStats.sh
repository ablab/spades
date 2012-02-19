#!/bin/sh

echo "The size of the clustered index is " `wc -l data/input/clustered.prd`
echo "The size of the etalon index is " `wc -l data/input/etalon.prd`
java -cp ./build GenStartStats
cd data/input
mv test_tp.prd tp.prd
mv test_fp.prd fp.prd

sort -rnk 4,4 fp.prd > fpr.prd
mv fpr.prd fp.prd
sort -rnk 4,4 tp.prd > tpr.prd
mv tpr.prd tp.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd

cd ../..
java -cp ./build GenStartStats > /dev/null

cd data/input
mv test_fp.prd fn.prd
rm test_tp.prd

sort -rnk 4,4 fn.prd > fnr.prd
mv fnr.prd fn.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd
cd ../..
