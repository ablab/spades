#!/usr/bin/env perl

$in = shift @ARGV;
$minGap =  shift @ARGV;
open(IN, $in) or die "cannot open $in\n";
$cur = <IN>;
$cur =~ /(^\d+) /;
$curPos = $1;
$maxDiff = 0;
while(<IN>) {
  $next = $_;
  $next =~ /(^\d+)/;
  $nextPos = $1;
  if ($nextPos - $curPos > $minGap) {
    print "$nextPos $curPos\n";
  }
  if ($nextPos - $curPos > $maxDiff) {
				$maxDiff = $nextPos - $curPos;
  }
  $curPos = $nextPos;
}
print "$maxDiff\n";
