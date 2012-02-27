#!/usr/bin/env perl

$in = shift @ARGV;

open(IN, "$in");


$lastone = 0;
while (<IN>) {
  $line = $_;
  $line =~ /\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(-?\d+)/;
  if ($1 eq "1") {
    $lastone = $line;
  }
}

print "lastone: $lastone\n";
