#!/usr/bin/env perl

if ($#ARGV < 0) {
		print "usage: ComputeCoverage.pl readMapFile \n";
		print "  -begin regionBegin start the genome at 'regionBegin' (lowest read)\n";
		print "  -end regionEnd end coverage at 'regionEnd' (higest read end)\n";
		print "  -cov coverage  Consider coverage to start and end at 'coverage'.\n";
		print "  -map mapFile   Write output to 'mapFile'\n";
	print 
  exit(0);
}

$readMapFile = shift @ARGV;

$coverage = 1;
$regionBegin = 0;
$regionEnd = -1;
$mapFile = "";

while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-cov") {
    $coverage = shift @ARGV;
  } elsif ($opt eq "-map") {
    $mapFile = shift @ARGV;
  }
	elsif ($opt eq "-begin") {
			$regionBegin = shift @ARGV;
	}
	elsif ($opt eq "-end" ){
			$regionEnd = shift @ARGV;
	}
}

open(RMF, $readMapFile) or die "cannot open $readMapFile\n";

@covered = ();
$begin = $regionBegin;
$end   = $regionEnd;
for ($i = $begin; $i < $end; $i++ ) {
  $covered[$i - $begin] = 0;
}

$totLength = 0;
$nMapped   = 0;
$maxReadEnd= 0;
while(<RMF>) {
		if ($_ =~ /(\d+) (\d+)/) {
				$start = $1; $end = $2;
				
				if ((($start > $regionStart &&
							$start < $regionEnd) ||
						 ($end > $regionStart &&
							$end < $regionEnd) ||
						 ($start < $regionStart &&
							$end > $regionEnd)) ||
						$regionEnd == -1){
						$nMapped++;
						#find read startpoints
						if ($start < $regionStart) {
								$covStart = $regionStart;
						}
						else {
								$covStart = $start;
						}
						# find read endpoints
						if ($regionEnd != -1 &&
								$end > $regionEnd) {
								$covEnd = $regionEnd;
						}
						else {
								$covEnd = $end;
						}
						if ($covEnd > $maxReadEnd) {
								$maxReadEnd = $covEnd;
						}
#						print "covering $covStart ... $covEnd\n";
						for ($c = $covStart; $c <= $covEnd; $c++ ) {
								$covered[$c - $regionStart]++;
						}
						$totLength += $covEnd - $covStart + 1;
				}
		}
}

# now compute the number of contigs
$nContigs = 0;
$contigStart = -1;
if ($covered[0] != 0) {
  $nContigs++;
  $contigStart = 0;
}

$totCovered = 0;
$totMult;
@zeroStretch = ();
@contigLengths = ();
@contigStarts = ();
@contigEnds   = ();
$lastZero     = 0;
if ($covered[0] >= $coverage) {
		push @contigStarts, 0;
}

# If no read end was specified on the command line, guess that a read covers
# the end of the sequence.
if ($regionEnd == -1) {
		$regionEnd = $maxReadEnd;
}

for ($p = 1; $p <= $regionEnd - $regionStart; $p++ ) {
  if ($covered[$p] < $coverage)  {
    # not considered a contig. If this is adjacent to a spot
    # that was a contig, then we have a break in coverage
    if ( $covered[$p-1] >= $coverage) {
      push @zeroStretch, 1;
      $lastZero = $#zeroStretch;
      push @contigLengths, $p - $contigStart;
      push @contigEnds, $p;
    }
    else {
      $zeroStretch[$lastZero]++;
    }
  }
  # Maybe this is the start of a contig?
  if ($covered[$p] >= $coverage &&
      $covered[$p-1] < $coverage) {
    $contigStart = $p;
    push @contigStarts, $p;
    $nContigs++;
  }
  if ($covered[$p] >= $coverage) {
    $totCovered++;
    $totMult += $covered[$p];
  }
}

$avgMult = $totMult / $totCovered;
$avgLength = $totLength / $nMapped;
print "total covered: $totCovered  mapped reads: $nMapped\n";
print "$nContigs with average coverage: $avgMult avg mapped length: $avgLength $nMapped mapped\n";
print "stretches of length 0 are: @zeroStretch \n";

print "The contig coords are: \n";
$nLt500 = 0;
$largeCov = 0;
for ($i = 0; $i < $#contigStarts; $i++ ) {
		$len = $contigEnds[$i] - $contigStarts[$i];
  print "$contigStarts[$i] ... $contigEnds[$i] ($len)\n";
		if ($len < 500) {
				$nLt500++;
		}
		else {
				$largeCov += $contigEnds[$i] - $contigStarts[$i] + 1;
		}
}
$tot = scalar @contigStarts;
print "$nLt500 less than 500 (out of $tot)\n";
print "total coverage is $totCovered, large contigs (>500) cover: $largeCov\n";

if ($mapFile ne "") {
		open(MAP, ">$mapFile") or die "cannot write to $mapFile\n";
		$last = $#covered;
		for ($i = 0; $i <= $last; $i++) {
				print MAP "$covered[$i]\n";
		}
}
