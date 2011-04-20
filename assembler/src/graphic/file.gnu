#! /usr/bin/gnuplot -persist

set term png enhanced font "/usr/share/fonts/truetype/msttcorefonts/couri.ttf,8" \
xffffff x000000 size 640,480

set encoding iso_8859_1    		
set tmargin 1        			
set tics out        			
set size 1.0,1.0			
set nokey				
set grid				
set xzeroaxis lt -1			
set yzeroaxis
 
plot "./daily.dat" using 1:2 smooth csplines linestyle 1
