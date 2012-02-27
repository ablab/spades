#!/usr/bin/env perl

if ($#ARGV < 0) {
		print "usage: ComputeOverlapCoverage.pl readMapFile \n";
		print "  -begin 'begin'   Start the genome at 'regionBegin' (lowest read)\n";
		print "  -end 'end'     Stop computing coverage at 'regionEnd' (higest read end)\n";
		print "  -vertexSize v  Find breaks in overlaps for vertex size 'v'\n";
		print "  -cov coverage  Consider coverage to start and end at 'coverage'.\n";
		print "  -map mapFile   Write output to 'mapFile'\n";
		print "  -repeatMap map Use coordinates in 'map' to define hard-stops on \n";
		print "                 repeat coverage.\n";
	print 
  exit(0);
}

$readMapFile = shift @ARGV;

$coverage = 1;
$regionBegin = 0;
$regionEnd = -1;
$mapFile = "";
$repeatMapFile = "";
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
	elsif ($opt eq "-vertexSize" ) {
			$vertexSize = shift @ARGV;
	}
	elsif ($opt eq "-repeatMap") {
			$repeatMapFile = shift @ARGV;
	}
	else {
			print "bad option: $opt\n";
			exit(1);
	}
}

open(RMF, $readMapFile) or die "cannot open $readMapFile\n";

while(<RMF>) {
		$_ =~ /(\d+) (\d+)/;
		$start = $1;  $end = $2;
		push @readMap, [($start, $end)];
}

if ($repeatMapFile ne "") {
		open(RPMF, $repeatMapFile) or die "cannot open $repeatMapFile\n";
		while(<RPMF>) {
				$_ =~ /(\d+) (\d+)/;
				$start = $1;  $end = $2;
				push @repeatMap, [($start, $end)];
		}
}

# now sort this
print "sorting ";
@sortedMap = sort TwoDSort @readMap;
print "done\n";


$curCoverageBegin = $sortedMap[0][0];
$curCoverageEnd   = $sortedMap[0][1];


$maxReadEnd =  0;
$prevCovered = 0;
@readsInCurContig = ();
@curContigLength  = ();
$curContig = 0;
for ($i = 1; $i < scalar @sortedMap; $i++) {
		$prevCovered = 0;
		
		if ($sortedMap[$i-1][1] - $sortedMap[$i][0] > $vertexSize) {
#				print "ovp: $sortedMap[$i-1][0] $sortedMap[$i-1][1] $sortedMap[$i][0] $sortedMap[$i][1]\n";
				# There is an overlap for this read as well, mark it.
				for ($pos = $sortedMap[$i-1][0]; $pos < $sortedMap[$i-1][1]; $pos++) {
						$covered[$pos]++;
				}
				$prevCovered = 1;
				@readsInCurContig[$curContig]++;
		}
		else {
				# There is an overlap for this read as well, mark it.
				if ($prevCovered) {
						for ($pos = $sortedMap[$i-1][0]; $pos < $sortedMap[$i-1][1]; $pos++) {
								$covered[$pos]++;
						}
						@readsInCurContig[$curContig]++;
						$curCoverageEnd = $osrtedMap[$i-1][1];
				}
				@curContigLength[$curContig] = $curContigEnd - $curContigBegin + 1;
				$curContig++;
				$curContigBegin = $sortedMap[$i][0];

		}
						
		if ($maxReadEnd < $sortedMap[$i][1]) {
				$maxReadEnd = $sortedMap[$i][1];
		}
}

# Force breaks at the boundaries of repeats
for ($i = 0; $i < scalar @repeatMap; $i++ ){
		for ($j = $repeatMap[$i][0]; $j < $repeatMap[$i][1]; $j++) {
				$covered[$j] = 0;
				$covered[$j] = 0;
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

$regionStart = $sortedMap[0][0];
# If no read end was specified on the command line, guess that a read covers
# the end of the sequence.
if ($regionEnd == -1) {
		$regionEnd = $maxReadEnd;
}

print "checking form 1 to $regionEnd\n";
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

$avgMult   = $totMult   / $totCovered;
print "total covered: $totCovered.\n";
$nContigLengths = scalar @contigLengths;
print "$nContigs ($nContigLengths) with average coverage: $avgMult avg mapped length: $avgLength $nMapped mapped\n";
print "stretches of length 0 are: @zeroStretch \n";

print "The contig coords are: \n";
$nLt500 = 0;
$nGt500 = 0;
$largeCov = 0;
for ($i = 0; $i < $#contigStarts; $i++ ) {
		$len = $contigEnds[$i] - $contigStarts[$i];
		print "$contigStarts[$i] ... $contigEnds[$i] $contigLengths[$i] $readsInCurContig[$i]\n";
#		push @contigLengths, $len;
		if ($len < 500) {
				$nLt500++;
		}
		else {
				$largeCov += $contigEnds[$i] - $contigStarts[$i] + 1;
				$nGt500++;
		}
}

# calc N50 size
@css = sort {$a <=>$b} @contigLengths;

$totSize = 0;
for ($i = 0; $i < @css; $i++) {
	$totSize = $totSize + @css[$i];
}
$cumSize = 0;
$n50Index = -1;
for ($i = $#css; $i >= 0 && $n50Index == -1; $i--) {
  $cumSize = $cumSize + @css[$i];
  if ($cumSize > ($totSize / 2)) {
     $n50Index = $i;
  }
}

$n50 = $css[$n50Index];

$tot = scalar @contigStarts;
print "$nLt500 less than 500 (out of $tot)\n";
print "Total coverage is $totCovered, $nGt500 large contigs (>500) cover: $largeCov\n";
if ($n50 != -1) {
		print "the N50 size is $n50\n";
}


sub TwoDSort {
#		print "in mycmp, @{$a}[0],  @{$b}[0]\n";
		@{$a}[0] <=> @{$b}[0];
}
