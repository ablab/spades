#!/usr/bin/env perl

$fFile = shift @ARGV;
$rFile = shift @ARGV;

open(F, $fFile) or die "cannot open $fFile\n";
open(R, $rFile) or die "cannot open $rFile\n";

$nF = 0;
while (<F>) {
		chomp;
		if (/^>(.*)_F3/) {
				$fTitle = $1;
				$fTitles{$fTitle} = $nF;
				$nF++;
		}
}

$nR = 0;
while (<R>) {
		chomp;
		if (/^>(.*)_R3/) {
				$rTitle = $1;
				$rIndex = $nF + $nR;
				if (exists $fTitles{$rTitle}) {
						print "$rTitle". "_F3 " . $rTitle . "_R3 $nF $rIndex 0\n";
				}
				$nR++;
		}
}
