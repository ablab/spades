#!/usr/bin/env perl
if ($#ARGV < 3) {
  print "usage: ComputeErrorProfile.pl minErrors seqLength file1.aln file2.aln...\n";
  exit(1);
}
$minErrors = shift @ARGV;
$seqLen = shift @ARGV;
$nSeq = 0;
@errCount = ();
$totReads = ();
for ($p = 0; $p < $seqLen; $p++) {push @errCount, 0;}
		
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
				if ($nerr == $minErrors || $minErrors == 0) {
						if ($strand == 1) {
								$pre = $alnStr;
								$alnStr = reverse($alnStr);
						}
						@aln = split(//, $alnStr);
						for ($a = 0; $a <= $#aln; $a++) {
								$errCount[$a] += $aln[$a];
						}
						$nSeq++;
				}
		}
}

for ($e = 0; $e <= $#errCount; $e++) {
		$rate = $errCount[$e] / $nSeq;
		print "$rate ";
}
