#! /bin/tcsh -f
if ($#argv < 3) then
  echo "usage: $0 infile outfile port [minscore]"
  exit
endif
if ($#argv > 3) then
  echo gfClient out=blast8 -minScore=$4 darth.ucsd.edu $3 /  $1 $2 
  /home/mchaisso/bin/alpha/gfClient -out=blast8  -minScore=$4 darth.ucsd.edu $3 / $1 $2 
else
  /home/mchaisso/bin/alpha/gfClient darth.ucsd.edu $3 / $1 $2 -out=blast8 
endif
 
