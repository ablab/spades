#!/bin/bash

for file in $(ls $1/saves/)
do 
	echo $file
	cmp  $1/saves/$(basename $file) $2/saves/$(basename $file)
done

