#!/usr/bin/env perl

$exeDir = "~/projects/mcsrc";
if ($#ARGV < 4 ) {
  print "usage: PrintReadsInInterval.pl reads.map start end out\n";
  print "  reads.map is a file of mapped read positions created by first \n";
  print "  blasting all reads to a finished sequence, then using the script\n";
  print "  assembly/MapReads.pl to print a map of read_title start end\n";
  exit(1);
}

$readsName    = shift @ARGV;
$readPosName = shift @ARGV;
$start = shift @ARGV;
$end   = shift @ARGV;
$outFile = shift @ARGV;


system("$exeDir/UTILS/PrintRange.pl $readPosName $start $end > $readsName.$start.$end");
system("$exeDir/UTILS/ExtractSequences.pl $readsName.$start.$end $readsName > $outFile");
      
