#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 locations_file output_file
$locationsFile = shift @ARGV;

open(LOCSIN, "$locationsFile") or die "Cannot open $locationsFile\n";
