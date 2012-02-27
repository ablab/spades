#!/usr/bin/env csh
echo "$1 $2"
$1 $2
set outfile = `First.pl $2`
set s = $status
@ nFailed = 0
if (! -e $outfile) then
  echo "failed command " >> failed_jobs.txt
  echo "$1 $2" >> failed_jobs.txt
  echo "REASON: no output" >> failed_jobs.txt
  @ nFailed = 1
  # try it just once more
  $1 $2
endif
if (-z $outfile) then
  echo "failed command " >> failed_jobs.txt
  echo "$1 $2" >> failed_jobs.txt
  echo "REASON: zero output size" >> failed_jobs.txt
  # try it just once more
  $1 $2
  @ nFailed = $nFailed + 1
endif
set d = `date`
echo "$d $1 $2 STATUS: $status NTRIES: $nFailed" >> command_log.txt
  
