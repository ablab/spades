#!/usr/bin/env perl

while (<>) {
  $_ =~ /(\S+) (\d+) (\d+)/;
  #print $_;
  $region = $1;
  $start  = $2;
  $end    = $3;
  $d = $start - $end;
  #print "$start, $end, $d \n";
  if (abs($d) > 50000) {
    print "$d - $region $start $end \n";
  }
}
