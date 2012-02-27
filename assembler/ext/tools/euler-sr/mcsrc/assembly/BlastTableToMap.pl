#!/usr/bin/env perl

if ($#ARGV < 2) {
  print "usage: BlastHitsToMap.pl blastTableFile readsFile mapFile\n";
  exit(0);
}
$blastTableFile = shift @ARGV;
$readsFile    = shift @ARGV;
$mapFile      = shift @ARGV;

`RetainHighestBlastHit.pl $blastTableFile > $blastTableFile.best`;
`sort $blastTableFile.best > $blastTableFile.best.sorted`;
`grep ">" $readsFile > $readsFile.names`;
`sort $readsFile.names > $readsFile.names.sorted`;
`MapReads.pl $readsFile.names.sorted $blastTableFile.best.sorted`;
