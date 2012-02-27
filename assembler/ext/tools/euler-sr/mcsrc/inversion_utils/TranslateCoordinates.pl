#!/usr/bin/env perl

$in = shift @ARGV;

open (IN, $in) or die "cannot open $in\n";
while (<IN>) {
  $line = $_;
  $line =~ /(.*):(\d+)\-(\d+)/;
  $chr = $1; $start = $2; $end = $3;

}

