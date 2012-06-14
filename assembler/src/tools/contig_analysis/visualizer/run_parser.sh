#!/bin/bash

for file in $(ls *.stdout)
do
	python parser.py $file >> list.txt
done
