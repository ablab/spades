#!/usr/bin/env perl
if ($#ARGV < 2) {
    print "Computes coverage of the genome using positions stored in\n";
    print "reads.  Any coverage that drops below 'thresh' initiates\n"; 
    print "a gap in coverage.\n\n";
		print "usage: ComputeCoverage.pl reads thresh vertexSize\n";
		exit(0);
}
$in = shift @ARGV;
$thresh = shift @ARGV;
$v  = shift @ARGV;
$printLengths = 0;
$readLength = -1;
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq  "-l") {
				$printLengths = 1;
		}
		elsif ($opt eq "-rl") {
				$readLength = shift @ARGV;
		}
}

@genome = ();
open(IN, "$in") or die "cannot open $in\n";
while(<IN>) {
		if ($_ =~ /^>(\d+)_(...)_(\d+)/) {
				$dir = $2;
				$pos = $3;
		}
		elsif ($_ =~ /pos=(\d+).*strand=(\d+)/) {
				$dir = $2;
				$pos = $1;
		}
		else {
				if ($readLength == -1) {
						$rl = length ($_);
				}
				else {
						$rl = $readLength;
				}
				if ($dir eq "REV") {
						for ($i = 0; $i < $rl - $v + 1; $i++) {
								$genome[$pos + $rl - $i - $v]++;
						}
				}
				else {
						for ($i = 0; $i < $rl - $v + 1; $i++) {
								$genome[$pos + $i]++;
						}
				}
		}
}

$stop = 0;
$start = 0;
for ($p = 1; $p < $#genome ; $p++) {
		if ($genome[$p] >= $thresh && $genome[$p+1] < $thresh) {
				$stop = $p;
				if ($printLengths == 0) {
						print "$start $stop\n";
				}
				else {
						$len = $stop - $start;
						print "$len\n";
				}
		}
		elsif ($genome[$p] >= $thresh && $genome[$p-1] < $thresh) {
				$start = $p ;
		}
}
