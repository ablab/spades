#!/bin/bash

for i in 21 33 55
do
echo "Generating k.hpp for K="$i
./gen_k.sh $i 

echo "Compiling assembler for K="$i
make rd

echo "Running assembler for K="$i
run rd > "K$i_iteration.log"
done
