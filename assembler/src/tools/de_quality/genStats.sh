#!/bin/sh
javac GenStartStats.java
java GenStartStats
mv test_tp.prd tp.prd
mv test_fp.prd fp.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd

java GenStartStats
mv test_fp.prd fn.prd
rm test_tp.prd

mv etalon.prd temp.prd
mv clustered.prd etalon.prd 
mv temp.prd clustered.prd
