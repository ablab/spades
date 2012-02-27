#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: $0 coordsList geneList\n";
		exit(0);
}

$coordsListFile = shift @ARGV;
$geneListFile = shift @ARGV;

open(CF, "$coordsListFile") or die "cannot open $coordsListFile\n";

@contigStart = ();
@contigEnds  = ();
while(<CF>) {
		$_ =~ /(\d+) (\d+) (\d+)/;
		$start = $1;
		$end   = $2;
		push @contigStart, $1;
		push @contigEnd, $2;
}

open(GL, "$geneListFile") or die "cannot open $geneListFile\n";
# the first 3 lines are descriptions
<GL>;
<GL>;
<GL>;
$nGenes = 0;
$nContained = 0;
while(<GL>) {
		$_ =~ /(\d+)..(\d+)/;
		$start = $1;
		$end   = $2;
		$contained = 0;
		$c = 0;
		while ($c <= $#contigStart && $contained == 0) {
				if ($start >= $contigStart[$c] &&
						$end <= $contigEnd[$c]) {
						$contained = 1;
						$nContained++;
				}
				$c++;
		}
		$nGenes++;
}
				
print "$nContained / $nGenes\n";
