#!/bin/bash 

cd ../../../
./gen_k 55
make rdt 
./run rdt 
cd src/test/teamcity/

