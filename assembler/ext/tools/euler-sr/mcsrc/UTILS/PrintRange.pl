#!/usr/bin/env perl
if ($#ARGV < 2) {
  print "Print read titles within a range\n";
  print "usage: PrintRange.pl reads.map start end [-coords]\n";
  exit 1;
}
 
$readPosIn = shift @ARGV;
$start = shift @ARGV;
$end   = shift @ARGV;
$printCoords = 0;
while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-coords") {
    $printCoords = 1;
  }
}

open(RPI, $readPosIn) or die "cannot open $readPosIn\n";

while(<RPI>) {
  $line = $_;
  $line =~ /(\S+) (\d+) (\d+)/;
  $title = $1;
  $seqStart = $2;
  $seqEnd   = $3;
  if (($seqStart >= $start && $seqStart <= $end) ||
      ($seqEnd  >= $start && $seqEnd <= $end) ||
      ($seqStart <= $start && $seqEnd >= $end)) {
    print "$title";
    if ($printCoords) {
      print " $seqStart $seqEnd";
    }
    print "\n";
  }
}
