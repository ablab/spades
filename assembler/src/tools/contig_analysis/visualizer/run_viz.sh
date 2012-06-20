#!/bin/bash

rm list.txt
touch list.txt

for file in $(ls *.stdout)
do
	python parser.py $file >> list.txt
done

python misassemblies.py -o plot -a list.txt --arcs

rm list.txt
