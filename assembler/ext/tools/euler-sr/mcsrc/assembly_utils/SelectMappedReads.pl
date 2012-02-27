#!/usr/bin/env perl
if ($#ARGV < 2) {
		print "usage: SelectMappedReads.pl read.map lowBound highBOund\n";
		exit(1);
}

$in = shift @ARGV;
$lowBound = shift @ARGV;
$highBound = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";

while(<IN>) {
		$_ =~ /(\d+) (\d+) (\S+) (\d)/;
		$start = $1; $end = $2;
		if (($start >= $lowBound &&
				 $end <= $highBound)) {
				print ">$3\n";
		}
}
