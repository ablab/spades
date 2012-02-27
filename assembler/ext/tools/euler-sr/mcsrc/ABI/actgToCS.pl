#!/usr/bin/env perl

$in = shift @ARGV;
$prim = shift @ARGV;
open(IN, $in) or die "cannot open $in\n";
while(<IN>) {
  if (/^#/) {
  }
  elsif (/^>/) {
    print;
  }
  else {
    /.(.*)/;
    $end = $1;
    chomp($end);
    $end =~ tr/ACTG/0123/;
    print "$prim$end\n";
  }
}
