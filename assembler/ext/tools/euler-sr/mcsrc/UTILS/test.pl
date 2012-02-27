#!/usr/bin/env perl
use POSIX;
# spawn a child process.

$st = "hello joe";
$sub1 = substr($st,0,1);
$sub2 = substr($st,-1,1);
$sub3 = substr($st, 0, -2);
print "got strings '$sub1', '$sub2', '$sub3'\n";



sub getLoadAvg {
  my ($machine) = @_;
  $topstr = `ssh $machine  -n -t -q \"top -n 1 \"`;
  @topstr = split /\n/, $topstr;
  foreach $s (@topstr) {
    if ($s =~ /load averages:  (\d+\.\d+)/) { 
      return $1;
    }
    elsif ($s =~ /load average: (\d+\.\d+)/) {
      return $1;
    }
  }
  return -1;
}
  
