#!/usr/bin/env perl
use POSIX;

if ($#ARGV < 1) {
		print "usage: FindGapsInCoverage.pl in kmer [minLength]\n";
		exit(0);
}

$in   = shift @ARGV;
$kmer = shift @ARGV;
$minLength = 0;
if ($#ARGV >= 0) {
		$minLength = shift @ARGV;
}
open(IN, "$in") or die "cannot open $in\n";

$curPos = 0;
$curLine = <IN>;

@coverage = ();

$readIndex = 0;
while(<IN>) {
		$curLine = $_;
		$curLine =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d)/;
		$curStart = $1; $curEnd = $2; $curTitle = $3; $curStrand = $4;
		
		$readLength = $curEnd - $curStart + 1;
		for($p = $curStart; $p < $curEnd - $kmer; $p++) {
				@coverage[$p]++;
		}
		++$readIndex;
#		if ($readIndex % 10000 == 9999) {
#				print ".";
#		}
#		if ($readIndex % 500000 == 499999) {
#				print "\n";
#		}
}
print "\n";

# now compute the number of contigs

$lenCovered = scalar @coverage;
$inContig = 1;
@contigStarts = ();
@contigEnds   = ();
$pos = 0;
while($pos < $lenCovered && @coverage[$pos] <= 1) {$pos++;}

while($pos < $lenCovered) {
		# scan to end of contig
		push @contigStarts, $pos;
		while($pos < $lenCovered && @coverage[$pos] > 1) {$pos++;}
		push @contigEnds, $pos;
		
		# scan to beginning of next contig
		while ($pos < $lenCovered && @coverage[$pos] <= 1){ $pos++;}
}
			
# print the contigs
$nContigs = scalar @contigEnds;
if ($nContigs > 0) {
		$covered = 0;
		$meanLength = 0;
		for ($c = 0; $c < $#contigStarts; $c++) {
				$covered += $contigEnds[$c] - $contigStarts[$c] + 1;
				$contigLength = $contigEnds[$c] - $contigStarts[$c] + 1;
				if ($contigLength > $minLength) {
						print "$contigLength\n";
				}
#				push @contigLengths, $contigEnds[$c] - $contigStarts[$c] + 1;
		}
		$contigLenghsLine = join(" ", @contigLengths);
		$meanLength = $covered / $nContigs;

#		print "total coverage: $covered in $nContigs average: $meanLength\n";
#		print "N50: $n50\n";
}
