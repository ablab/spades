#!/usr/bin/env csh
if ($# < 3) then
  echo "Usage: RunBlastz.csh ref qry \"args\" output"
  exit 1
endif

blastz $1 $2 $4 > $3

if ($status) then
  date >> ~/failed_blastz.txt
  echo "blastz $1 $2 $4 > $3" >> ~/failed_blastz.txt
endif

@ failed = 0
if (-z $3) then
  date >> ~/failed_jobs.txt
  pwd >> ~/failed_jobs.txt
  echo "blastz $1 $2 > $3" >> ~/failed_jobs.txt
  echo "PROBLEM: zero output size" >> ~/failed_jobs.txt
 @ failed = 1
endif

if ($failed != 1) then
  tail -1 $3 | grep -q "eof"
  if ($status) then
    date >> ~/failed_jobs.txt
    pwd >> ~/failed_jobs.txt
    echo "blastz $1 $2 $4 > $3" >> ~/failed_jobs.txt
    echo "PROBLEM: truncated output" >> ~/failed_jobs.txt
    @ failed = 1
  endif
endif

# report an error if any occurred

@ numtries = 0
while ($failed == 1 && $numtries < 3) 
  blastz $1 $2 $4 > $3
  @ failed = 0
   if (-z $3) then
     date >> ~/failed_jobs.txt
     pwd >> ~/failed_jobs.txt
     echo "blastz $1 $2 $4 > $3" >> ~/failed_jobs.txt
     echo "PROBLEM: zero output size" >> ~/failed_jobs.txt
     @ failed = 1
   endif
  if ($failed != 1) then
    tail -1 $3 | grep -q "eof"
    if ($status) then
      date >> ~/failed_jobs.txt
      pwd >> ~/failed_jobs.txt
      echo "blastz $1 $2 $4> $3" >> ~/failed_jobs.txt
      echo "PROBLEM: truncated output" >> ~/failed_jobs.txt
      @ failed = 1
    endif
  endif
  @ numtries = $numtries + 1
end
echo "running: blastz $1 $2 $4 > $3  STATUS: $failed  TRIES: $numtries" >> blastz_report.txt
  
  # try to re-run the script
  
