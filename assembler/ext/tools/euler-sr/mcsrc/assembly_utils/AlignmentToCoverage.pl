#!/usr/bin/env perl

if ($#ARGV < 1) {
		print "usage: $0 alignmentFile minVotes\n";
		exit(0);
}

$in       = shift @ARGV;
$wordSize = shift @ARGV;


open(IN, $in) or die "cannot open $in\n";

while(<IN>) {
		$line = $_;
		$line =~ /(\d+)\s+(\d+)\s+([01]*)\s+(\d+).*/;
		$pos = $0;
		$alnStr = $3;
		@aln = split(//, $alnStr);
		$end = $#aln;
		for ($p = 0; $p <= $end; $p++) {
		
