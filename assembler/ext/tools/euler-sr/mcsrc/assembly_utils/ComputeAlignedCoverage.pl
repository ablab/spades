#!/usr/bin/env perl
if ($#ARGV < 3) {
  print "usage: ComputeErrorProfile.pl minErrors seqLength file1.aln file2.aln...\n";
  exit(1);
}
$nSeq = 0;
@errCount = ();
$totReads = ();

@genome = ();
$nReads = 0;
foreach $file (@ARGV) {
		open(IN, $file) or die "cannot open $file\n";
		while(<IN>) {
				#parse lines like this:
				#30089   1       00000101000100000000000000000000000000000000000000 3
				$_=~/(\d+)\s+(\d+)\s+(\d+)\s+(\d+).*/;
				$pos = $1;
				$strand = $2;
				$alnStr = $3;
				$nerr = $4;
				@genome[$pos]++;
				$nReads++;
		}
}

$sumsq = 0;
$window = 8;
for($p = 0; $p <= $#genome - $window; $p++) {
		$winSum = 0;
		for ($w = $p; $w < $p + $window; $w++) {
				$winSum += $genome[$w];
		}
		$sumsq += ($winSum * $winSum);
}
$mean = ($nReads/$#genome) * $window;
$var = $sumsq/($#genome) - $mean*$mean ;
$genomeLen = $#genome;
print "#reads: $nReads genomeLen: $genomeLen mean: $mean var: $var\n";

