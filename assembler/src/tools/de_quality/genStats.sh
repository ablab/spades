#!/bin/sh
java -cp ./build GenStartStats
cd data/input
mv test_tp.prd tp.prd
mv test_fp.prd fp.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd

cd ../..
java -cp ./build GenStartStats > /dev/null

cd data/input
mv test_fp.prd fn.prd
rm test_tp.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd
cd ../..
