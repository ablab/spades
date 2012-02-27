#!/usr/bin/env csh
if ($#argv < 1) then
  echo "Usage: PrintStats.csh  edgeFile "
  exit(1)
endif
set baseFile = $1:r
set N50500 =  `pcl $1  -minLength 500 | N50`
echo "N50 500     $N50500";
set N50 = `pcl $1 | N50`
echo "N50         $N50";
set med = `pcl $1 | median`
echo "median      $med";
setenv len `pcl $1 -minLength 500 | add.pl`
@ longlen = $len / 2
echo "len (>500)  $longlen"
setenv len `pcl $1 | add.pl`
@ len = $len / 2;
echo "len (all)   $len"
set nEdge = `pcl $1 | wc -l`
echo "#edges      $nEdge"
set nLongEdge = `pcl $1 -minLength 500 | wc -l`
echo "#long edges $nLongEdge"

set numSources = `${EUSRC}/assembly/${MACHTYPE}/printGraphSummary $baseFile -sources | wc -l`
echo "# sources   $numSources"


