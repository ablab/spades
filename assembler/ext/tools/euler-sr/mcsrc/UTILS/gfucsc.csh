#! /bin/tcsh -f
if ($#argv > 2) then
  echo gfClient darth.ucsd.edu 9998 /  $1 $2 -out=blast8 -minScore=$3
  gfClient darth.ucsd.edu 9998 / $1 $2 -out=blast8  -minScore=$3
else
  gfClient darth.ucsd.edu 9998 / $1 $2 -out=blast8 
endif
