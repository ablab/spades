#!/usr/bin/env perl
if ($#ARGV < 0) {
		print "usage: PrintUniqueEdgeIntervalCoords.pl intvFile (minEdgeLength)\n";
}
$intvFile = shift @ARGV;
open(IN, $intvFile) or die "cannot open $intvFile\n";
$minEdgeLength = 0;
if ($#ARGV >= 0) {
		$minEdgeLength = shift @ARGV;
}

while(<IN>) {
		if ($_ =~ /EDGE.*Multiplicity 1/) {
				$intv = <IN>;
#parse:
#				INTV 0 734064 40 0
				$intv =~ /INTV\s+(\d+)\s+(\d+)\s+(\d+).*/;
				$strand = $1;
				$pos    = $2;
				$len    = $3;
				$start = $pos;
				$end   = $pos + $len;
				if ($strand == 0 && $len >= $minEdgeLength) {
						print "$start $end $len\n";
				}
		}
}
				
