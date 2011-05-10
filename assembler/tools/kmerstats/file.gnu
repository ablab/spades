#!/usr/bin/gnuplot -persist

set terminal png
set output "outfile.png"
plot 'daily.dat' using 1:2 notitle with lines
