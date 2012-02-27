#!/usr/bin/env csh
if ($#argv == 0) then
  echo "usage: uunpc.csh file";
  exit 1
endif

perl -pi -e "s/\r//" $1
