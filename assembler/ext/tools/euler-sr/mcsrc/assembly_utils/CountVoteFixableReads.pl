#!/usr/bin/env perl

if ($#ARGV < 1) {
		print "usage: $0 alignmentFile minVotes\n";
		exit(0);
}

$in = shift @ARGV;
$minVotes = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
# parse lines like this:
#17735   1       00000000000000000000100000000010000000000010000000 3    GATAAAAAGGACTTTGAGGGCTCAGGGGAGGAGGGGTAGGGAGGTGAGAG

$totalFixable = 0;
$total = 0;
$allOk = 0;
$middleError = 0;
while(<IN>) {
		$line = $_;
		$line =~ /(\d+)\s+(\d+)\s+([01]*)\s+(\d+).*/;
		$alnStr = $3;
		@aln = split(//, $alnStr);
		$end = $#aln;
		$fixable = 1;
		if ($aln[0] == 1 or $aln[$end] == 1) {
#				print "$aln[0] $aln[$end]\n";
				$fixable = 0;
		}
		else {
				$onePos = -1;
				$p = 0;
				$thisOk = 1;
				$nerr = 0;
				for ($p = 0; $p <= $end; $p++) {
						if ($aln[$p] == 1) {
								if ($p - $onePos - 1 < $minVotes) {
										# found two nearly adjacent errors, can't
										# fix this 
										$fixable = 0;
										$thisOk = 0;
								}
								$onePos = $p;
								$nerr++;
						}
				}
#				$allOk += $thisOk;
				if ($nerr == 0) {
						$allOk++;
				}
				# check to see if there is one error
				if ($nerr == 1) {
						$mi = 0;
						for ($p = 0; $p <= $end && $mi == 0; $p++) {
								if ($aln[$p] == 1 && $p > $minVotes && $p < ($end - $minVotes)) {
										$middleError++;
										$mi = 1;
								}
						}
				}
		}
		$totalFixable += $fixable;
		$total++;
}

print "$totalFixable / $total, $allOk perfect $middleError\n";
