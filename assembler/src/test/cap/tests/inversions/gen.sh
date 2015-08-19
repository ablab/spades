#!/bin/bash

numtest=0

function gentest() {
  ((numtest++));
  echo "$1 $2 tests/${numtest}f.fasta tests/${numtest}r.fasta" | ./gen > tests/$numtest;
}

echo "Creating tests";

rm gen &> /dev/null;
g++ -O2 -Wall -Wextra -o gen gen.cpp
if [ $? -ne 0 ]; then
  exit 1;
fi

mkdir tests &> /dev/null;

gentest 1000 1
gentest 1000 2
gentest 50000 50

rm gen
echo $numtest > tests/summary
