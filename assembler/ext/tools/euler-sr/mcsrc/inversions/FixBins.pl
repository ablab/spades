#!/usr/bin/env perl

if ($#ARGV != 0) {
  print "usage: FixBins.pl binfile\n";
  print "removes 'query:' from the binfile";
  exit(0);
}
$infile = shift @ARGV;
