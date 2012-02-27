#!/usr/bin/env perl
use strict;
use POSIX;

if ($#ARGV != 3) {
  printf("usage: $0 inmatefile readlocfile outmatefile readvar\n");
  exit(-1);
}

my($inmate, $readloc, $outmate, $pctvar) = @ARGV;

open(INMATE, "$inmate") or die "cannot open $inmate\n";
open(READLOC, "$readloc") or die "cannot open $inmate\n";
open(OUTMATE, ">$outmate") or die "cannot open $inmate\n";


my @locations = ();
my $temp;
while (<READLOC>) {
  push @locations, split(" ", $_);
}
close(READLOC);
my($fstart, $fend, $rstart, $rend, $dist, $var);

while(<INMATE>) {
  /([\w\.]+) (\d+) (\d+)/;
  $fstart = shift(@locations);
  $fend   = shift(@locations);
  $rstart = shift(@locations);
  $rend   = shift(@locations);
  $dist   = $rend - $fstart + 1;
  $var    = POSIX::floor($dist * $pctvar / 2);
  printf(OUTMATE "$1 $2 $3 %d %d\n", $dist - $var, $dist + $var);
}
close(INMATE);
close(OUTMATE);


  
