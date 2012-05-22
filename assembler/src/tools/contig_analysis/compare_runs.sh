#!/bin/bash

for file in $(ls $1/saves/*.sqn)
do
	echo $(basename $file)
	python compare_fasta.py $1/saves/$(basename $file)  $2/saves/$(basename $file)
done
