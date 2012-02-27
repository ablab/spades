#!/usr/bin/env perl
if ($#ARGV != 1) {
  print "usage: $0 infasta outbins\n";
  print "prints the coordinates as parsed from the fasta file\n";
  exit(0);
}

$in = shift @ARGV;
$out = shift @ARGV;

open (IN, $in) or die "cannot open $in\n";
open (OUT, ">$out") or die "cannot open $out\n";

print OUT "bin 0\n";
while (<IN>) {
  $line = $_;
  if ($line =~ /(\w+)\.(\d+)\.(\d+)/) {
    $spec = $1; $start = $2; $end = $3;
    print OUT "$spec $start $end\n";
  }
}
