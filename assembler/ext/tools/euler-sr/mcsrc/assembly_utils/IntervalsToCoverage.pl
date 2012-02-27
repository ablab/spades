#!/usr/bin/env perl

$intvFile = shift @ARGV;
$covFile = shift @ARGV;

open(INTV, $intvFile) or die "cannot open $intvFile\n";
open(COV, ">$covFile") or die "cannot open $covFile\n";

@cov = ();
while(<INTV>) {
		if (/EDGE (\d+) Length (\d+) Multiplicity (\d+)/) {
				if ($#cov > 0) {
						print COV "#edge $edge length: $length\n";
						for ($c = 0; $c <= $#cov; $c++) {
								print COV "$c $cov[$c]\n";
						}
						print COV "\n\n";
				}
				@cov = ();
				$edge = $1;
				$len  = $2;
				$mult = $3;
				for ($i = 0; $i < $len; $i++ ){
						$cov[$i] = 0;
				}
		}
		else {
				$_ =~ /INTV (\d+) (\d+) (\d+) (\d+)/;
				$read = $1;
				$readPos = $2;
				$readLen = $3;
				$edgePos = $4;
				$edgePosEnd = $edgePos + $readLen;
				for ($e = $edgePos; $e < $edgePosEnd; $e++) {
						$cov[$e]++;
				}
		}
}
				
