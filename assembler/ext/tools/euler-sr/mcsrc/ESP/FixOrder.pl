#!/usr/bin/env perl

$in = shift @ARGV;

open(IN, "$in") or die "cannot open $in\n";

@coords = ();
while(<IN>) {
  @vals = split(/\s+/, $_);
  if ($vals[0] < $vals[1]) {
    $first = $vals[0];
    $second = $vals[1];
  }
  else {
    $first = $vals[1];
    $second = $vals[2];
  }
  print  "$first $second $vals[2] $vals[3] $vals[4]\n";
}
