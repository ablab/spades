#!/usr/bin/env perl
if ($#ARGV < 1) {
  print "usage: TrimReadEnds.pl in out\n";
  exit(0);
}
$in = shift @ARGV;
$out = shift @ARGV;

open(IN, "$in") or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";

while(<IN>) {
if ($_~ /^>/) {
  $title = $_;
  $readSeq = "";
}
else {
  $line = $_;
  chomp $line;
  $readSeq .= $line;
  }
}
