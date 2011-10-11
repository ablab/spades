#!/bin/bash


function iterative_run
{
echo "Iterative run started"


for i in 21 33 55
do

echo "Generating k.hpp for K="$i
./gen_k.sh $i 

echo "Compiling assembler for K="$i
make rd

echo "Running assembler for K="$i
./run rd > $1/"K$i.log"

done

echo "Iterative run finished"
}

if [ $# = 0 ]
then 
    set $1 "."
fi

if [ $1 = "." ]
then
    echo "Warning: writing logs to the root directory"
else
    echo "Creating $1 log directory"
    mkdir $1
fi

iterative_run $1 > $1/overall.log
