#!/usr/bin/env perl
use POSIX;
if ($#ARGV < 1) {
		print "usage: divfasta.pl infile num_files [out_base]\n";
		exit(1);
}

$inFile = shift @ARGV;
$numFiles = shift @ARGV;
$inFile =~ /(.*)\.[^\.]*$/;
$outBase = $1;
if ($#ARGV >= 0) {
		$outBase = shift @ARGV;
}
# count the number of entries first.
open(IN, $inFile) or die "cannot open $inFile\n";
$numSeq = 0;
while(<IN>) {
		if (/^>/)  {
				$numSeq++;
		}
}

$curSeq = 0;
close(IN);
open(IN, $inFile) or die "cannot open $inFile\n";
$numSeqPerFile = POSIX::floor($numSeq / $numFiles);
$outNumber = 0;
$outFile = $outBase . ".$outNumber.fasta";
open(OUT, ">$outFile") or die "cannot open $outFile\n";
$seq = "";
while(<IN>) {
		if (/^>/)  {
				if ($seq ne "") {
						print OUT "$seq";
				}
				$curSeq++;
				if ($curSeq > $numSeqPerFile and $outNumber < $numFiles-1) {
#       print OUT $seq;
						close OUT;
						$outNumber++;
						$outFile = $outBase . ".$outNumber.fasta";
						open(OUT, ">$outFile") or die "cannot open $outFile\n";
						$curSeq = 0;
				}
				$seq = $_;
		}
		else {
				$seq .= $_;
		}
}

print OUT "$seq";

