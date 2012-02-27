#!/usr/bin/env perl
if ($#ARGV < 3) {
		print "usage: SelectReadsByCoverage.pl readsFile genome vertexSize minCoverage\n";
		exit(1);
}

$readsFile = shift @ARGV;
$genome    = shift @ARGV;
$vertexSize = shift @ARGV;
$minCoverage= shift @ARGV;
$genomeSize = `pcl $genome`;
chomp $genomeSize;
$totalUncovered = 0;
@genomeCov = ();
for ($g = 0; $g < $genomeSize; $g++) {
		push @genomeCov, 0;
		$totalUncovered++;
}

$maxGaps = 50;

open(RF, $readsFile) or die "cannot open $readsFile\n";


# Don't sweat some of the boundaries.

$bound = 50;
for ($p = 0; $p < $bound; $p++) {
		$genomeCov[$p] = $minCoverage;
		$totalUncovered--;
}
for ($p = $#genomeCov; $p >= $#genomeCov - $bound; $p--) {
		$genomeCov[$p] = $minCoverage;
		$totalUncovered--;
}
$readIndex = 0;
while(<RF>) {
		if ($_ =~ /^>(\d+)_(...)_(\d+)/) {
				$dir = $2;
				$pos = $3;
				$title = $_;
				$read = <RF>;
				$len = length($read)-1;
#				print "$pos $len $vertexSize\n";
				for ($p = $pos; $p < $pos + $len - $vertexSize + 1; $p++) {
						$genomeCov[$p]++;
						if ($genomeCov[$p] == $minCoverage) {
								$totalUncovered--;
						}
				}
				print $title;
				print $read;
		}
		if ($totalUncovered == 0) {
				exit(0);
		}
		if ($readIndex % 1000 == 0) {
				print STDERR "$readIndex $totalUncovered\n";
				$numNotCovered = 0;
				for ($tp = 0; $tp <= $#genomeCov && $numNotCovered < $maxGaps; $tp++) {
						if ($genomeCov[$tp] < $minCoverage) {
								$numNotCovered++;
						}
				}
				if ($numNotCovered < $maxGaps) {
						exit(0);
				}
		}
		$readIndex++;
}

