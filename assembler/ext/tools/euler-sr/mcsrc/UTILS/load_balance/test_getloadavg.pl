#!/usr/bin/env perl

$loadAvg   = getLoadAvg();
print "$loadAvg\n";

sub getLoadAvg {
  $os = `uname`;
  chomp $os;
  if ($os eq "FreeBSD") {
    $la = `getloadavg.bsd`;
    print "la: $la\n";
  }
  else {
    $la = `getloadavg`;
  }
  $la =~ /([\d\.]+)\s+.*/;
  return $1;
}
