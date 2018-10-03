#!/bin/sh

total="$(grep ">" $1 | wc -l)"
mapped="$(grep ">" $2 | wc -l)"
echo $mapped
echo $total
echo $((mapped*100/total))
